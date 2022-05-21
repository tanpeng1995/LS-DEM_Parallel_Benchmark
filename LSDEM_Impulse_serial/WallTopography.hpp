#ifndef WALLTOPOGRAPHY_HPP_
#define WALLTOPOGRAPHY_HPP_
#include "definitions.hpp"
#include "Grain3d.hpp"
#include "utilities.hpp"

class WallTopography{
public:
  WallTopography(){}
  ~WallTopography(){}
  WallTopography(const string & topoFile, const double & angle, const Vector3d & origin, const double & mu, const double & kn){
    _rotMatrix << cos(angle), 0., -sin(angle),
                      0.,     1.,     0.,
                  sin(angle), 0., cos(angle);
    _origin = origin; // _origin << 0.,0.,600.;
    _mu     = mu;
    _kn     = kn;
    _ks     = 0.9 * _kn;
    ifstream file(topoFile.c_str());
  	if(file.fail()){
      cout << "ERROR: Topography file does not exist." << endl;
    }
  	string  line;
  	string 	partial;
  	istringstream iss;
    /* dimension: x: 2000, y: 800, z: 150 */
    getline(file, line);
  	iss.str(line);
  	getline(iss, partial, ' ');
  	_xdim = atoi(partial.c_str());
  	getline(iss, partial, ' ');
  	_ydim = atoi(partial.c_str());
  	getline(iss, partial, ' ');
  	_zdim = atoi(partial.c_str());
  	iss.clear();
    cout<<"dimX: "<<_xdim<<" dimY: "<<_ydim<<" dimZ: "<<_zdim<<endl;
    /* level set */
    while (getline(file, line)){
      _lset.push_back(atof(line.c_str()));
    }
    if(_xdim * _ydim * _zdim != _lset.size()){
      cout << "ERROR: levelset dimension does not match. Actually get: "<<_lset.size()<< endl;
    }
  	file.close();
    cout<<"global rotation matrix of topography is: "<<endl;
    cout<<_rotMatrix<<endl;
  }
  bool bCircleCheck(const std::unique_ptr<Grain3d> & grain) const{
    Vector3d point = _rotMatrix * ( grain->getPosition() - _origin );
    return grain->getRadius() > getLSetValue(point);
  }
  bool bTempCircleCheck(const std::unique_ptr<Grain3d> & grain) const{
    Vector3d point = _rotMatrix * ( grain->getTempPosition() - _origin );
    Vector3d pos = grain->getTempPosition();
    Vector3d disp = grain->getTempPosition() - _origin;
    return grain->getRadius() > getLSetValue(point);
  }
  void findGrainWallContact(std::unique_ptr<Grain3d> & grain, double dt){
    Vector3d normalForce(0.,0.,0.);
    Vector3d shearForce(0.,0.,0.);
    Vector3d grainMoment(0.,0.,0.);
    Vector3d normal;
    Vector3d point;
    double   penetration = 0.;
    double   minPenetration = std::numeric_limits<double>::max();
    Vector3d minNormal(0.,0.,0.);
    Vector3d ptcm;
    Vector3d v, ds;
    vector<Vector3d> pointList = grain->getPointList();
    for (int i = 0; i < pointList.size(); ++i) {
      Vector3d grainForce(0.,0.,0.);
      point = _rotMatrix * ( pointList[i] - _origin );
      findPenetration(point, penetration, normal);
      if (penetration < 0) {
        normal = _rotMatrix.transpose() * normal;
        ptcm = pointList[i] - grain->getTempPosition() - penetration * normal;
        normalForce = -penetration*normal*_kn;
        v = grain->getVelocity() + grain->getOmegaGlobal().cross(ptcm);
        ds = (v-v.dot(normal)*normal)*dt;
        shearForce = -ds*_kn;
        if(shearForce.norm() > _mu*normalForce.norm()){
          double forceMagnitude = _mu*normalForce.norm();
          shearForce = shearForce/shearForce.norm() * forceMagnitude;
        }
        grainForce = normalForce + shearForce;
        grainMoment += ptcm.cross(grainForce);
        if(penetration < minPenetration){
          minPenetration = penetration;
          minNormal = normal;
        }
      }
      /* end penetration < 0. */
    }
    /* change velocity */
    if(minPenetration < 0.){
      grain->addVelocity(-minPenetration / dt * minNormal);
    }
    /* change omega */
    if(grainMoment.norm() > 0.){
      Vector3d dOmega(0.,0.,0.);
      Vector3d omega = grain->getOmega();
      Vector3d omegaN = omega;
      Vector3d momentInertia = grain->getMomentInertia();
      double   gDamping = 0.5;
      Vector3d principleMoment = grain->getTempRotMatrix().transpose()*grainMoment;
      for (int i = 0; i < 3; i++) {
    		dOmega(0) = (principleMoment(0) + omegaN(1)*omegaN(2)*(momentInertia(1)-momentInertia(2)) - gDamping*momentInertia(0)*omegaN(0) )*dt/momentInertia(0);
    		dOmega(1) = (principleMoment(1) + omegaN(2)*omegaN(0)*(momentInertia(2)-momentInertia(0)) - gDamping*momentInertia(1)*omegaN(1) )*dt/momentInertia(1);
    		dOmega(2) = (principleMoment(2) + omegaN(0)*omegaN(1)*(momentInertia(0)-momentInertia(1)) - gDamping*momentInertia(2)*omegaN(2) )*dt/momentInertia(2);
    		omegaN = omega + dOmega/2.;
    	}
      //if( dOmega.norm() > 1. ) dOmega = dOmega / dOmega.norm();
      /* this is important at the first step, initial overlaps, may not need it, if the initial positions are good */
      omega += dOmega;
      grain->changeOmega(omega);
    }
    /* end findGrainWallContact */
  }

  void findGrainWallContactImpulse(std::unique_ptr<Grain3d> & grain, double dt, vector<ContactInfo> & contactList){
    Vector3d normalForce(0.,0.,0.);
    Vector3d shearForce(0.,0.,0.);
    Vector3d grainMoment(0.,0.,0.);
    Vector3d normal;
    Vector3d point;
    Vector3d ptcm;
    Vector3d v, ds;
    Vector3d netNormal(0.,0.,0.), netBV(0.,0.,0.), netV(0.,0.,0.);
    int      nContact = 0;
    double   penetration = 0.;
    vector<Vector3d> pointList = grain->getPointList();
    for (int i = 0; i < pointList.size(); ++i) {
      Vector3d grainForce(0.,0.,0.);
      point = _rotMatrix * ( pointList[i] - _origin );
      findPenetration(point, penetration, normal);
      if (penetration < 0) {
        normal = _rotMatrix.transpose() * normal;
        ptcm = pointList[i] - grain->getPosition() - penetration * normal;
        v = grain->getVelocity() + grain->getOmegaGlobal().cross(ptcm);
        netNormal += normal; netBV += ptcm; netV += v; nContact++;
      }
      /* end penetration < 0. */
    }
    if(nContact > 0){
      netNormal /= nContact; netBV /= nContact; netV /= nContact;
      double relVel = netNormal.dot(netV);
      if(relVel < -threshold){
        contactList.emplace_back(grain->getId(), -1, netBV, Vector3d(0.,0.,0.), netNormal, netV, relVel, 0.);
      }
    }
  }

  private:
    double getLSetValue(const Vector3d & point) const{
      double x = point(0); double y = point(1); double z = point(2);
    	if (x+1 > (double)_xdim || y+1 > (double)_ydim || z+1 > (double)_zdim || x < 0 || y < 0 || z < 0 ){
        return std::numeric_limits<double>::max();
      }
    	return getGridValue(round(x),round(y),round(z));
    }
    double getGridValue(const int & x, const int & y, const int & z) const {
    	return (double)_lset[z*_ydim*_xdim + y*_xdim + x];
    }
    bool findPenetration(const Vector3d & point, double & penetration, Vector3d & normal) const {
      double x = point(0); double y = point(1); double z = point(2);
    	if (x+1 > (double)_xdim || y+1 > (double)_ydim || z+1 > (double)_zdim || x < 0 || y < 0 || z < 0 )
    		return false;
    	// check if the point is close to the surface, if not, return false
    	if (getGridValue(round(x),round(y),round(z)) > 1.)
    		return false;
    	// if the point is close to the surface, find the value by performing trilinear interpolation
    	int x0 	= floor(x);
    	int y0 	= floor(y);
    	int z0 	= floor(z);
    	int x1 	= ceil(x);
    	int y1 	= ceil(y);
    	int z1 	= ceil(z);
    	double p000 = getGridValue(x0,y0,z0);
    	double p001 = getGridValue(x0,y0,z1);
    	double p010 = getGridValue(x0,y1,z0);
    	double p011 = getGridValue(x0,y1,z1);
    	double p101 = getGridValue(x1,y0,z1);

    	double xm 	= getGridValue(x1,y0,z0) - p000;
    	double ym 	= p010 - p000;
    	double zm	  = p001 - p000;
    	double xym	= -xm - p010 + getGridValue(x1,y1,z0);
    	double xzm	= -xm - p001 + p101;
    	double yzm	= -ym - p001 + p011;
    	double xyzm = -xym + p001 - p101 - p011 + getGridValue(x1,y1,z1);
    	double dx 	= x - double(x0);
    	double dy 	= y - double(y0);
    	double dz	  = z - double(z0);
    	penetration = p000 + xm*dx + ym*dy + zm*dz + xym*dx*dy + xzm*dx*dz + yzm*dy*dz + xyzm*dx*dy*dz;
    	// if the point lies outside the surface, return false
    	if (penetration >= 0)
    		return false;
    	// finally, if the point lies in the surface, find the gradient, normalize, and return true
    	normal << xm + xym*dy + xzm*dz + xyzm*dy*dz,
    					    ym + xym*dx + yzm*dz + xyzm*dx*dz,
    					    zm + xzm*dx + yzm*dy + xyzm*dx*dy;
    	normal /= normal.norm();
    	return true;
    }

    Matrix3d _rotMatrix;
    Vector3d _origin;
    double   _mu;
    double   _kn;
    double   _ks;
    int      _xdim;
    int      _ydim;
    int      _zdim;
    vector<float> _lset;
};

static std::unique_ptr<WallTopography> wallTopography;

#endif /*WALLTOPOGRAPHY_HPP_*/
