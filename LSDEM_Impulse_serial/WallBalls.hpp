#ifndef WALLBALLS_HPP_
#define WALLBALLS_HPP_
#include "definitions.hpp"
#include "Grain3d.hpp"

class WallBalls {
public:
	// default constructor
	WallBalls() {
		_kn = 0;
		_kmn = 0.;
    _nvert = 0;
    _nangular = 0;
    _pressure = 0;
    _height = 0;
    _l0 = 0;
    _ballRadius = 0;
    _ballMass = 0;
    _segHeight = 0;
    _ballMoin = 0;
    _nballs = 0;
	}

	/* WallBalls for an inclined plane */
	// constructor if you actually want a working wallBalls
	WallBalls(const double & length, const double & height, const double & ballRadius, const double & density, const double & kn, const double & kmn, const double & angle, const Vector3d & pos, const int & startIdx){
		_pressure = 0.;
		_kn = kn;
		_kmn = kmn;
		_startIdx = startIdx;
		// how to arithmetic?
		double n = length/ballRadius/2.;
		_nangular = floor(n);
		double a = 1-cos(2.*M_PI/double(_nangular));
		_ballRadius = 1.1*ballRadius;
		_ballMass = 4./3.*M_PI*_ballRadius*_ballRadius*_ballRadius*density;
		_ballMoin = 0.4*_ballMass*_ballRadius*_ballRadius;
		_l0 = 2*_ballRadius;
		double htheta = M_PI/double(_nangular);
		// note: for some reason this doesn't quite lead to all springs being the same starting length fixme
		_segHeight = sqrt(3.) * _ballRadius;
    _nvert = ceil(height/_segHeight);
		_nballs = _nvert*_nangular;
		// set up vectors and zero out vectors that need to be zeroed
		_velocities.resize(_nballs);
		_positions.resize(_nballs);
		_ballForces.resize(_nballs);
		_shears.resize(_nballs);
		_normals.resize(_nballs);
		_height = (double)(_nvert-1)*_segHeight+2*_ballRadius;
		for (int i = 0; i < _nballs; i++) {
			_velocities[i] << 0., 0., 0.;
			_positions[i] << 0., 0., 0.;
			_ballForces[i] << 0., 0., 0.;
			for (int j = 0; j < 6; j++) {
				_shears[i][j] << 0., 0., 0.;
				_normals[i][j] << 0., 0., 0.;
			}
		}
		// fill ball vectors
		double sinTheta = sin(angle);
		double cosTheta = cos(angle);
		double x, y;
		int ball_id;
		for (int t = 0; t < _nangular; t++) {
			x = (double(t)/double(_nangular))*length;
			for (int h = 0; h < _nvert; h+=2) {
				ball_id = t*_nvert + h;
				y = _ballRadius+_segHeight*double(h);
				_positions[ball_id] << y*sinTheta, x, y*cosTheta;
				_positions[ball_id] += pos;
			}
			// offset theta by a half-angle for odd-numbered layers
			x = ((double(t)+.5)/double(_nangular))*length;
			for (int h = 1; h < _nvert; h+=2) {
				ball_id = t*_nvert + h;
				y = _ballRadius+_segHeight*double(h);
				_positions[ball_id] << y*sinTheta, x, y*cosTheta;
				_positions[ball_id] += pos;
			}
		}

		// fill faces vector
		_faces.resize( (_nvert-1)*2*(_nangular-1) );
		int upsidedown = 0; // are the faces upside-down?
		int t = 0; int h = 0; //t: horizontal h: height
		for (int f = 0; f < _faces.size(); f++) {
			if (upsidedown == 0 && h%2 == 0) { _faces[f] << t*_nvert+h, (t+1)*_nvert+h, t*_nvert+h+1; }
			else if (upsidedown == 1 && h%2 == 1){ _faces[f] << t*_nvert+h, (t+1)*_nvert+h-1, (t+1)*_nvert+h; }
			else if (upsidedown == 0 && h%2 == 1) { _faces[f] << t*_nvert+h, (t+1)*_nvert+h, (t+1)*_nvert+h+1; }
			else { _faces[f] << t*_nvert+h, t*_nvert+h-1, (t+1)*_nvert+h; } // upsidedown == 1 && h%2 == 0
			t = t == _nangular-2 ? 0 : t+1;
			//t = (t+1)%_nangular;
			if (t == 0) {
				if (upsidedown == 0) { h++;}
				upsidedown = 1-upsidedown;
			}
		}
		// fill neighbors vector
		_neighbors.resize(_nballs);
		for (int t = 1; t < _nangular-1; t++) {
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

		{
			int t = 0;
			for (int h = 0; h < _nvert; h++) {
				int n = t*_nvert+h;
				int nr = ((t+1)%_nangular)*_nvert+h;
				int nur = (n+1 + (h%2)*_nvert)%(_nvert*_nangular);
				int nul = (n+1 - ((h+1)%2)*_nvert + _nvert*_nangular)%(_nvert*_nangular);
				int nl = ((t+_nangular-1)%_nangular)*_nvert+h;
				int nbl = (n-1 - ((h+1)%2)*_nvert + _nvert*_nangular)%(_nvert*_nangular);
				int nbr = (n-1+(h%2)*_nvert)%(_nvert*_nangular);
				if (h==0) { _neighbors[n] << nr, nur, -1, -1, -1, -1; }
				else if (h==_nvert-1) {	_neighbors[n] << nr, nbr, -1, -1, -1, -1; }
				else { _neighbors[n] << nr, nur, nbr, -1, -1, -1; }
			}
		}

		{
			int t = _nangular-1;
			for (int h = 0; h < _nvert; h++) {
				int n = t*_nvert+h;
				int nr = ((t+1)%_nangular)*_nvert+h;
				int nur = (n+1 + (h%2)*_nvert)%(_nvert*_nangular);
				int nul = (n+1 - ((h+1)%2)*_nvert + _nvert*_nangular)%(_nvert*_nangular);
				int nl = ((t+_nangular-1)%_nangular)*_nvert+h;
				int nbl = (n-1 - ((h+1)%2)*_nvert + _nvert*_nangular)%(_nvert*_nangular);
				int nbr = (n-1+(h%2)*_nvert)%(_nvert*_nangular);
				if (h==0) { _neighbors[n] << nl, nul, -1, -1, -1, -1; }
				else if (h==_nvert-1) {	_neighbors[n] << nl, nbl, -1, -1, -1, -1; }
				else { _neighbors[n] << nl, nul, nbl, -1, -1, -1; }
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
	}

	void findWallForceMoment(const std::unique_ptr<Grain3d> & grain, vector<ContactInfo> & contactList){
		Vector3d grainPos = grain->getPosition();
		Vector3d ball2Grain, ball2GrainDir, ball2Lset, normal, velocity;
		double penetration;
		double GamaN = -2*sqrt(_kn*_ballMass*grain->getMass()/(_ballMass+grain->getMass()))*log(0.6)/sqrt(M_PI*M_PI+log(0.6)*log(0.6));
		for(int ball_id = 0; ball_id < _nballs; ++ball_id){
			ball2Grain = _positions[ball_id] - grainPos;
			ball2GrainDir = ball2Grain / ball2Grain.norm();
			ball2Grain -= _ballRadius * ball2GrainDir;
			if(ball2Grain.squaredNorm() < grain->getRsq()){
				ball2Lset = grain->getRotMatrix().transpose()*ball2Grain + grain->getCmLset();
				if(grain->getLset().findPenetration(ball2Lset, penetration, normal)){
					normal = -grain->getRotMatrix()*normal;
					velocity = grain->getVelocity()+grain->getOmegaGlobal().cross(ball2Grain)-_velocities[ball_id];
					//_ballForces[ball_id] += penetration*normal*_kn - GamaN*velocity.dot(normal)*normal;
					if(normal.dot(velocity) < -threshold){
						contactList.emplace_back(grain->getId(), _startIdx+ball_id, ball2Grain, normal*ballRadius, normal, velocity, normal.dot(velocity), 0.);
					}
				}
			}
		}
	}

	void takeTimeStep(double dt){
		// fix bottom and top
		for (int i = 0; i < _nangular; i++) {
			//_velocities[i*_nvert + _nvert-1] << 0,0,0;
			_velocities[i*_nvert] << 0,0,0;
		}
		//fix sides
		for(int i = 0; i < _nvert; i++) {
			_velocities[i] << 0,0,0;
			_velocities[(_nangular-1)*_nvert+i] << 0,0,0;
		}
		for (int ball_id = 0; ball_id < _nballs; ++ball_id) {
			_positions[ball_id] += dt * _velocities[ball_id];
			_ballForces[ball_id] << 0., 0., 0.;
		}
	}

	void correctMembraneVelocity(double dt){
		/* neighbor force */
		Vector3d force(0.,0.,0.);
		Vector3d diff(0.,0.,0.);
		Vector3d normal(0.,0.,0.);
		double   Fmag, sint, cost, dist = 0.;
		Vector3d v(0.,0.,0.);
		Vector3d k, shear, ds;
		double GamaN = -2*sqrt(_kmn*_ballMass/2.)*log(0.6)/sqrt(M_PI*M_PI+log(0.6)*log(0.6));

		for (int ball_id = 0; ball_id < _nballs; ++ball_id) {
			for (int i = 0; i < 6; ++i) {
				int neighbor = _neighbors[ball_id][i];
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
				shear -= ds*_kmn;
				_shears[ball_id][i] = shear;
				_normals[ball_id][i] = normal;
				_ballForces[ball_id] += shear;
			}
		}

		// pressure forces
		for (int f = 0; f < _faces.size(); f++) {
			Vector3i face = _faces[f];
			Vector3d crossp = (_positions[face(2)]-_positions[face(0)]).cross(_positions[face(1)]-_positions[face(0)]);
			Vector3d ballForce = crossp*_pressure/6.;
			for (int i = 0; i < 3; i++) {
				_ballForces[face(i)] += ballForce;
			}
		}

		for (int ball_id = 0; ball_id < _nballs; ++ball_id) {
			_velocities[ball_id] += _ballForces[ball_id] / _ballMass * dt;
		}

		/* apply damping and zero out forces */
		for (int ball_id = 0; ball_id < _nballs; ++ball_id) {
			_velocities[ball_id] *= (1 - gDamping*dt);
			_ballForces[ball_id] << 0.,0.,0.;
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

	double getMass() const{
		return _ballMass;
	}

	Vector3d getVelocity(int id) const{
		return _velocities[id];
	}

	Vector3d getPosition(int id) const{
		return _positions[id];
	}

	double getBallRadius() const{
		return _ballRadius;
	}

	int getNumWalls() const{
		return _nballs;
	}

	void clearVelocities(){
		for(int i = 0; i < _nballs; ++i){
			_velocities[i] << 0., 0., 0.;
		}
	}

	double getMoin() const{
		return _ballMoin;
	}

	void addVelocity(int ball_id, Vector3d dv){
		_velocities[ball_id] += dv;
	}

	const vector<Vector3d> & getPositions() const {
		return _positions;
	}

	const vector<Vector3d> & getForces() const {
		return _ballForces;
	}

	const vector<Vector3d> & getVelocities() const {
		return _velocities;
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
	void changeStartIdx(const int & startIdx){
		_startIdx = startIdx;
	}
	const int getStartIdx() const{
		return _startIdx;
	}

private:
	Matrix3d qToR(const Vector4d & quat) {
		Matrix3d rotMatrix;
		rotMatrix(0,0) = -quat(0)*quat(0) + quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
		rotMatrix(0,1) = -2*(quat(0)*quat(1) - quat(2)*quat(3));
		rotMatrix(0,2) =  2*(quat(1)*quat(2) + quat(0)*quat(3));
		rotMatrix(1,0) = -2*(quat(0)*quat(1) + quat(2)*quat(3));
		rotMatrix(1,1) =  quat(0)*quat(0) - quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
		rotMatrix(1,2) = -2*(quat(0)*quat(2) - quat(1)*quat(3));
		rotMatrix(2,0) =  2*(quat(1)*quat(2) - quat(0)*quat(3));
		rotMatrix(2,1) = -2*(quat(0)*quat(2) + quat(1)*quat(3));
		rotMatrix(2,2) = -quat(0)*quat(0) - quat(1)*quat(1) + quat(2)*quat(2) + quat(3)*quat(3);
		return rotMatrix;
	}

	// must be updated together
	double _height;		// total height of wall (bottom is always at z=0 deal w/ it)
	double _segHeight;	// height of segment (starts at ball height but decreases during compression);

	// must be updated together
	vector<Vector3d> _velocities; // velocity of each ball
	vector<Vector3d> _positions;	// (x,y,z) position of each ball
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
	double _kmn;      // ball-ball stiffness
	double _ballRadius;	// radius of each ball
	double _ballMass;		// mass of each ball
	double _ballMoin;		// moment of inertia of each ball

	int _startIdx;
};

static std::vector<std::unique_ptr<WallBalls>> wallBalls(2);
static int num_balls = 0;

#endif /* WALLBALLS_HPP_ */
