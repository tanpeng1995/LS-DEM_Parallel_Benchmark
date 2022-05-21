/*
 * WALLBALLS_H_
 *
 *  Created on: Nov 17, 2015
 *      Author: Reid Kawamoto (Caltech)
 *  Modified on: June 21, 2020
 *      Author: Peng TAN (Berkeley)
 */

#ifndef WALLBALLS_H_
#define WALLBALLS_H_
#include "definitions.h"
#include "Grain3d.h"
#include "WallCap.h"

extern int openmpThreads;
extern int dynamicSchedule;

class WallBalls {
public:
	// default constructor
	WallBalls() {
		_kn = 0;
    _gdamping = 0;
    _nvert = 0;
    _nangular = 0;
    _pressure = 0;
    _height = 0;
    _l0=0;
    _kms = 0;
    _ballRadius=0;
    _ballMass=0;
    _segHeight=0;
    _kmn=0;
    _nballs=0;
	}

	// constructor if you actually want a working wallBalls
	WallBalls(const double & height, const double & initRadius, const double & ballRadius, const double & pressure,
    const double & density, const double & kn, const double & kmn, const double & kms, const double & gdamping):
		_pressure(pressure), _kn(kn), _kmn(kmn), _kms(kms), _radius(initRadius){
		_gdamping = gdamping;
		double n = M_PI*(initRadius+ballRadius)/ballRadius;
		_nangular = floor(n);
		double a = 1-cos(2.*M_PI/double(_nangular));
		_ballRadius = initRadius*(a+sqrt(2*a))/(2.-a);
		_ballMass = 4./3.*M_PI*_ballRadius*_ballRadius*_ballRadius*density;
		_l0 = 2*_ballRadius;
		double htheta = M_PI/double(_nangular);
		_segHeight = sqrt(4.*_ballRadius*_ballRadius-(initRadius+_ballRadius)*(initRadius+_ballRadius)*((1.-cos(htheta))*(1.-cos(htheta))+sin(htheta)*sin(htheta))); //WHY?
		_nvert = ceil(height/_segHeight);
		_nballs = _nvert*_nangular;
		_velocities.resize(_nballs);
		_positions.resize(_nballs);
		_ballForces.resize(_nballs);
		_touched.resize(_nballs);
		_shears.resize(_nballs);
		_normals.resize(_nballs);
		_rotMatrix << 1., 0., 0., 0., 1., 0., 0., 0., 1.;

		for (int i = 0; i < _nballs; i++) {
			_velocities[i] << 0., 0., 0.;
			_positions[i] << 0., 0., 0.;
			_ballForces[i] << 0., 0., 0.;
			_touched[i] = false;
			for (int j = 0; j < 6; j++) {
				_shears[i][j] << 0., 0., 0.;
				_normals[i][j] << 0., 0., 0.;
			}
		}

		for (int t = 0; t < _nangular; t++) {
			double curtheta = (double(t)/double(_nangular))*2.*M_PI;
			double x = (initRadius+_ballRadius)*cos(curtheta);
			double y = (initRadius+_ballRadius)*sin(curtheta);
			for (int h = 0; h < _nvert; h+=2) {
				int ball_id = t*_nvert + h;
				_positions[ball_id] << x, y, _ballRadius+_segHeight*double(h);
			}
			curtheta = ((double(t)+.5)/double(_nangular))*2.*M_PI;
			x = (initRadius+_ballRadius)*cos(curtheta); // change this to offset center
			y = (initRadius+_ballRadius)*sin(curtheta); // change this to offset center
			for (int h = 1; h < _nvert; h+=2) {
				int ball_id = t*_nvert + h;
				_positions[ball_id] << x, y, _ballRadius+_segHeight*double(h);
			}
		}

		_height = (double)(_nvert-1)*_segHeight+2*_ballRadius;
		moveHeight(height-_height);

		_faces.resize( (_nvert-1)*2*_nangular );
		int upsidedown = 0; // are the faces upside-down?
		int t = 0; int h = 0; //t: horizontal h: height
		for (int f = 0; f < _faces.size(); f++) {
			if (upsidedown == 0 && h%2 == 0) { _faces[f] << t*_nvert+h, ((t+1)%_nangular)*_nvert+h, t*_nvert+h+1; }
			else if (upsidedown == 1 && h%2 == 1){ _faces[f] << t*_nvert+h, ((t+1)%_nangular)*_nvert+h-1, ((t+1)%_nangular)*_nvert+h; }
			else if (upsidedown == 0 && h%2 == 1) { _faces[f] << ((t+0)%_nangular)*_nvert+h, ((t+1)%_nangular)*_nvert+h, ((t+1)%_nangular)*_nvert+h+1; }
			else { _faces[f] << t*_nvert+h,  ((t+0)%_nangular)*_nvert+h-1, ((t+1)%_nangular)*_nvert+h; } // upsidedown == 1 && h%2 == 0
			t = (t+1)%_nangular;
			if (t == 0) {
				if (upsidedown == 0) { h++;}
				upsidedown = 1-upsidedown;
			}
		}
		// fill neighbors vector
		_neighbors.resize(_nballs);
		for (int t = 0; t < _nangular; t++) {
			for (int h = 0; h < _nvert; h++) {
				int n = t*_nvert+h;
				int nr = ((t+1)%_nangular)*_nvert+h;
				int nur = (n+1 + (h%2)*_nvert)%(_nvert*_nangular);
				int nul = (n+1 - ((h+1)%2)*_nvert + _nvert*_nangular)%(_nvert*_nangular);
				int nl = ((t+_nangular-1)%_nangular)*_nvert+h;
				int nbl = (n-1 - ((h+1)%2)*_nvert + _nvert*_nangular)%(_nvert*_nangular);
				int nbr = (n-1+(h%2)*_nvert)%(_nvert*_nangular);
				if (h==0) { _neighbors[n] << nr, nur, nul, nl, -1, -1; }
				else if (h==_nvert-1) {	_neighbors[n] << nr, nl, nbl, nbr, -1, -1; }
				else { _neighbors[n] << nr, nur, nul, nl, nbl, nbr; }
			}
		}
		/* set normals */
		for (int ball_id = 0; ball_id < _nballs; ++ball_id) {
			for (int i = 0; i < 6; ++i) {
				int neighbor = _neighbors[ball_id](i);
				if (neighbor < 0) {break;}
				Vector3d diff = _positions[neighbor]-_positions[ball_id];
				_normals[ball_id][i] = diff/diff.norm();
			}
		}
	} // end constructor

	bool bCircleCheck(std::unique_ptr<Grain3d> & grain) const {
		Vector3d pos = grain->getPosition();
    if(sqrt(pos(0)*pos(0)+pos(1)*pos(1)) + grain->getRadius() < 0.8*_radius ) return false;
    return true;
  }

	void findWallForceMoment(std::unique_ptr<Grain3d> & grain, const double & dt, Vector6d & stressVoigt) {
		Vector3d grainForce(0.,0.,0.), grainMoment(0.,0.,0.);
		Vector3d grainPos = grain->getPosition();
		Vector3d branchVec, ball2Grain, ball2GrainDir, ball2Lset, normal, v, Fn;
		double penetration;
		double GamaN = -2*sqrt(_kn*_ballMass*grain->getMass()/(_ballMass+grain->getMass()))*log(0.6)/sqrt(M_PI*M_PI+log(0.6)*log(0.6));
		//GamaN = 0.;
		#pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads) private(branchVec, ball2Grain, ball2GrainDir, ball2Lset, normal, v, Fn, penetration) shared(grainPos, GamaN, _ballForces, stressVoigt) reduction(+:grainForce,grainMoment)
		for (int ball_id = 0; ball_id < _nballs; ++ball_id) {
			branchVec = grainPos - _positions[ball_id];
			ball2Grain = _positions[ball_id] - grainPos;
			ball2GrainDir = ball2Grain / ball2Grain.norm();
			ball2Grain -= _ballRadius * ball2GrainDir;
			if(ball2Grain.squaredNorm() < grain->getRsq()){
				ball2Lset = grain->getRotMatrix().transpose()*ball2Grain + grain->getCmLset();
				if(grain->getLset().findPenetration(ball2Lset, penetration, normal)){
					_touched[ball_id] = true;
					normal = grain->getRotMatrix()*normal;
					v = _velocities[ball_id]-grain->getVelocity()-grain->getOmegaGlobal().cross(ball2Grain);
					Fn = penetration*normal*_kn - GamaN*normal.dot(v)*normal;
					_ballForces[ball_id] -= Fn;
					grainForce  += Fn;
					grainMoment += ball2Grain.cross(Fn);
                                        /*
					stressVoigt(0) += Fn(0)*branchVec(0);
					stressVoigt(1) += Fn(1)*branchVec(1);
					stressVoigt(2) += Fn(2)*branchVec(2);
					stressVoigt(3) += 0.5*(Fn(1)*branchVec(2) + Fn(2)*branchVec(1));
					stressVoigt(4) += 0.5*(Fn(2)*branchVec(0) + Fn(0)*branchVec(2));
					stressVoigt(5) += 0.5*(Fn(1)*branchVec(0) + Fn(0)*branchVec(1));
                                        */
				}
			}
		}
		grain->addForce(grainForce);
		grain->addMoment(grainMoment);
	}

	void clearBallVelocities(){
		for(int i = 0; i < _nballs; ++i){
			_velocities[i] << 0.,0.,0.;
		}
	}

	double findVolume() const {
		double volume = 0;
		for (size_t h = 0; h < _nvert; h++) {
			for (size_t t = 0; t < _nangular; t++) {
				size_t idx0 = t*_nvert+h;
				size_t idx1 = ((t+1)%_nangular)*_nvert+h;
				volume += _positions[idx0](0)*_positions[idx1](1)-_positions[idx0](1)*_positions[idx1](0);
			}
		}
		return volume/double(_nvert)/2.*_height;
	}

	void takeTimeStep(double dt, string name){
		Vector3d force(0.,0.,0.);
		Vector3d diff(0.,0.,0.);
		Vector3d normal(0.,0.,0.);
		double   Fmag, sint, cost, dist = 0.;
		Vector3d v(0.,0.,0.);
		Vector3d k, shear, ds;
		double GamaN = -2*sqrt(_kmn*_ballMass/2.)*log(0.6)/sqrt(M_PI*M_PI+log(0.6)*log(0.6));
		//neighbor forces
		#pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads) private(force, diff, dist, v, normal, sint, cost, k, shear, ds, Fmag) shared(_neighbors, GamaN, _velocities)
		for (int ball_id = 0; ball_id < _nballs; ++ball_id) {
			for (int i = 0; i < 6; ++i) {
				int neighbor = _neighbors[ball_id](i);
				if (neighbor < 0) {break;}
				diff = _positions[neighbor]-_positions[ball_id];
				normal = diff/diff.norm();
				dist = diff.norm();
				v = _velocities[ball_id] - _velocities[neighbor];
				/* normal force */
				force = _kmn*(dist-_l0)*normal - GamaN*normal.dot(v)*normal;
				_ballForces[ball_id] += force;
				/* shear force */
				shear = _shears[ball_id][i];
				k = normal.cross(_normals[ball_id][i]);
				sint = k.norm(); cost = sqrt(1-sint*sint); if(isnan(cost)){ cost = 0.; }
				if (sint > 0){ k = k/(sint+std::numeric_limits<double>::epsilon()); }
				shear = shear*cost + k.cross(shear)*sint + k*k.dot(shear)*(1.-cost);
				/* update */
				ds = v * dt;
				shear -= ds*_kms;
				_shears[ball_id][i] = shear;
				_normals[ball_id][i] = normal;
				_ballForces[ball_id] += shear;
			}
		}
		//pressure force
		#pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads) shared(_faces, _positions, _pressure)
		for (int f = 0; f < _faces.size(); f++) {
			Vector3i face = _faces[f];
			Vector3d crossp = (_positions[face(2)]-_positions[face(0)]).cross(_positions[face(1)]-_positions[face(0)]);
			Vector3d ballForce = crossp*_pressure/6.;
			for (int i = 0; i < 3; i++) {
				_ballForces[face(i)] += ballForce;
			}
		}

		for (int i = 0; i < _nangular; i++) {
			// top layer
			int ball_id = i*_nvert + _nvert-1;
			_ballForces[ball_id] << 0,0,0;
			_velocities[ball_id] << 0,0,0;
			// bottom layer
			ball_id = i*_nvert;
			_ballForces[ball_id] << 0,0,0;
			_velocities[ball_id] << 0,0,0;
		}
    // time intergration
		for (int ball_id = 0; ball_id < _nballs; ++ball_id) {
			Vector3d vel = 1./(1.+_gdamping*dt/2.)*( (1.-_gdamping*dt/2.)*_velocities[ball_id] + dt*_ballForces[ball_id]/_ballMass );
			_velocities[ball_id] = vel;
			_positions[ball_id] += dt*_velocities[ball_id];
			_ballForces[ball_id] << 0., 0., 0.;
		}
	}

	void clearBallForcesMoments() {
		for (int i = 0; i < _ballForces.size(); i++) {
			_ballForces[i] << 0., 0., 0.;
		}
	}

	// finds the force at which the membrane pulls the top platen
	void computeInternalForceOnTopLayer(std::unique_ptr<WallCap> & wallCap, double dt) const {
		double GamaN = -2*sqrt(_kmn*_ballMass/2.)*log(0.6)/sqrt(M_PI*M_PI+log(0.6)*log(0.6));
		for (int i = 0; i < _nangular; i++) {
			int ball_id = i*_nvert + _nvert-1;
			Vector3d ramToBall = _positions[ball_id] - wallCap->getRamPos();
			Vector3d totalForce(0.,0.,0.);
			for (int n = 0; n < 6; n++) {
				int neighbor = _neighbors[ball_id][n];
				if (neighbor < 0) { break; }
				Vector3d diff = _positions[neighbor]-_positions[ball_id];
				double dist = diff.norm();
				Vector3d normal = diff / dist;
				Vector3d v = _velocities[ball_id] - _velocities[neighbor];
				Vector3d force = _kmn*(dist-_l0)*normal - GamaN*normal.dot(v)*normal;
				/* shear force */
				Vector3d shear = _shears[ball_id][n];
				Vector3d k = normal.cross(_normals[ball_id][n]);
				double sint = k.norm(); double cost = sqrt(1-sint*sint); if(isnan(cost)){ cost = 0.; }
				if (sint > 0){ k = k/(sint+std::numeric_limits<double>::epsilon()); }
				shear = shear*cost + k.cross(shear)*sint + k*k.dot(shear)*(1.-cost);
				Vector3d ds = v * dt;
				shear -= ds*_kms;
				totalForce = totalForce + force + shear;
			}
			wallCap->addForce(totalForce);
			wallCap->addMoment(ramToBall.cross(totalForce));
		}
	}

	// moves balls in top layer given by rotation R about pos
	void rotateTopLayer(const std::unique_ptr<WallCap> & wallCap) {
		Matrix3d R   = wallCap->getRotMatrix()*_rotMatrix.transpose();
		Vector3d pos = wallCap->getPosition();
		for (size_t t = 0; t < _nangular; t++) {
			size_t b = t*_nvert + _nvert-1;
			Vector3d trans = _positions[b] - pos;
			_positions[b] = R*trans/trans.norm()*_radius + pos;
		}
		_rotMatrix = wallCap->getRotMatrix();
	}

	void adjustHeights() {
		for (int h = 0; h < _nvert-1; h++) {
			double ringHeight = _segHeight*double(h) + _ballRadius;
			for (int t = 0; t < _nangular; t++) {
				int ball_id = t*_nvert+h;
				_positions[ball_id] << _positions[ball_id](0), _positions[ball_id](1), ringHeight;
			}
		}
	}


	// compresses wall vertically and moves balls correspondingly
	void moveHeight(const double & amount) {
		for (int i = 0; i < _nangular*_nvert; i++) {
			_positions[i](2) += amount*_positions[i](2)/_height;
		}
		_height += amount;
		_segHeight += amount/(double)_nvert;
	}


	void changeBallForces(int ball_id, Vector3d df){
		_ballForces[ball_id] += df;
	}

	void changel0(const double & l0) {
		_l0 = l0;
	}

	const double & getBallRadius() const {
		return _ballRadius;
	}
	const double & getKmn() const {
		return _kmn;
	}
	const double & getKms() const {
		return _kms;
	}
	const double & getL0() const {
		return _l0;
	}
	const double & getHeight() const{
		return _height;
	}
	const double & getSegHeight() const{
		return _segHeight;
	}
	const int & getNumWalls() const {
		return _nballs;
	}
	const int & getNangular() const {
		return _nangular;
	}
	const int & getNvert() {
		return _nvert;
	}
	const Vector3d & getPosition(const int & idx) const {
		return _positions[idx];
	}
	const vector<Vector3d> & getPositions() const {
		return _positions;
	}
	const vector<Vector3d> & getForces() const {
		return _ballForces;
	}
	const vector<Vector3i> & getFaces() const {
		return _faces;
	}
	const vector<Vector6i> & getNeighbors() const {
		return _neighbors;
	}
	const Vector3d & getVelocity(const int & idx) const {
		return _velocities[idx];
	}
	const vector<Vector3d> & getVelocities() const {
		return _velocities;
	}
	void changeVelocities(vector<Vector3d> ballVels){
		_velocities = ballVels;
	}
	void shrinkPressure(int n){
		_pressure = _pressure / n;
	}
	void changePosition(int bid, Vector3d pos){
		_positions[bid] = pos;
	}
	void changeVelocity(int bid, Vector3d vel){
		_velocities[bid] = vel;
	}
	void changeForce(int bid, Vector3d force){
		_ballForces[bid] = force;
	}
	vector<int> getBallsRank() {
		vector<int> ballsRank;
		for(int ball_id = 0; ball_id < _nballs; ++ball_id){
			if(_touched[ball_id]){
				ballsRank.push_back(ball_id);
				_touched[ball_id] = false;
			}
		}
		return ballsRank;
	}
private:
	// must be updated together
	double _height;		// total height of wall (bottom is always at z=0 deal w/ it)
	double _segHeight;	// height of segment (starts at ball height but decreases during compression);

	// must be updated together
	vector<Vector3d> _velocities; // velocity of each ball
	vector<Vector3d> _positions;	// (x,y,z) position of each ball
	vector<bool> _touched;
	vector<array<Vector3d, 6>> _shears;
	vector<array<Vector3d, 6>> _normals;

	// mesh-like stuff
	vector<Vector3d> _ballForces; // list of forces on the balls
	vector<Vector3i> _faces;		// list of faces (each face has the index of 3 balls that make up its vertices)
	vector<Vector6i> _neighbors;	// list of neighbors of each ball
	double _l0;							// initial spring length

	// these values are never changed
	int _nvert;			// number of balls in the vertical direction
	int _nangular;		// number of balls around
	int _nballs;		// total number of balls
	double _pressure;		// radial pressure (applied to each face)
	double _kn;				// ball normal stiffness
	double _kmn;			// mesh normal stiffness
	double _kms;			// mesh shear stiffness
	double _gdamping;		// global damping
	double _ballRadius;	// radius of each ball
	double _radius;
	double _ballMass;		// mass of each ball
	Matrix3d _rotMatrix;
};

#endif /* WALLBALLS_H_ */
