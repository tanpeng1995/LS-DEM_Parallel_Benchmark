/*
 * Grain3d.h
 *
 *  Created on: July 15, 2014
 *      Author: Reid Kawamoto(Caltech)
 *  Modified on: May 6, 2020
 *      Author: Peng TAN (Berkeley)
 *  Modified on: June 20, 2020
 *      Author: Peng TAN (Berkeley)
 */

#ifndef GRAIN3D_H_
#define GRAIN3D_H_

static int resolution;

class Grain3d{
public:
  Grain3d(){
    _radius   = 0;
    _kn       = 0;
    _ks       = 0;
    _cresN    = 0;
    _mass     = 0;
    _mu       = 0;
    _id       = 0;
    _density  = 0;
    _position << 0,0,0;
    _velocity << 0,0,0;
    _force    << 0,0,0;
    std::map<int, std::tuple<Vector3d, Vector3d>>().swap(_nodeFriction);
  }

  Grain3d(const Vector3d & position, const double & radius, const double & density, const double & kn, const double & ks, const double & mu, const int & gid){
    _mass = pow(0.95*radius,3.)*4./3.*M_PI;
    _force << 0.,0.,0.;
    _density = 1.;
    changeDensity(density);
    _radius   = 0.95*radius;
    _position = position;
    _id       = gid;
    _kn       = kn;
    _ks       = ks;
    _mu       = mu;
    _cresN    = 0.6;
    std::map<int, std::tuple<Vector3d, Vector3d>>().swap(_nodeFriction);
  }

  ~Grain3d(){
    _radius   = 0;
    _kn       = 0;
    _ks       = 0;
    _cresN    = 0;
    _mass     = 0;
    _mu       = 0;
    _id       = 0;
    _density  = 0;
    _position << 0,0,0;
    _velocity << 0,0,0;
    _force    << 0,0,0;
    std::map<int, std::tuple<Vector3d, Vector3d>>().swap(_nodeFriction);  }

  Grain3d(const Grain3d& other){
    _radius   = other._radius;
    _kn       = other._kn;
    _ks       = other._ks;
    _cresN    = other._cresN;
    _mass     = other._mass;
    _mu       = other._mu;
    _id       = other._id;
    _density  = other._density;
    _position = other._position;
    _velocity = other._velocity;
    _force    = other._force;
    _nodeFriction = other._nodeFriction;
  }

  Grain3d & operator=(const Grain3d& other){
    _radius   = other._radius;
    _kn       = other._kn;
    _ks       = other._ks;
    _cresN    = other._cresN;
    _mass     = other._mass;
    _mu       = other._mu;
    _id       = other._id;
    _density  = other._density;
    _position = other._position;
    _velocity = other._velocity;
    _force    = other._force;
    _nodeFriction = other._nodeFriction;
    return *this;
  }

  bool bCircleCheck(const Grain3d & other) const{
    Vector3d d = other.getPosition()-_position;
    double rsum = other.getRadius()+_radius;
    return d.norm() < rsum;
  }

  void findInterparticleForceMoment(Grain3d & other, const double & dt, Vector6d & stressVoigt){
    Vector3d force(0.,0.,0.);
    Vector3d newNodeNormal(0.,0.,0.);
    Vector3d newNodeShear(0.,0.,0.);
    int      newNodeContact;
    double   penetration = 0.;   // penetration amount
    Vector3d normal(0.,0.,0.);        // surface normal pointing out of other in the reference config, point toward *this
    Vector3d Fn(0.,0.,0.);            // normal force form a single point
    Vector3d v(0.,0.,0.);             // velocity of this grain wrt to other grain
    Vector3d Fs(0.,0.,0.);            // vector of frictional/shear force
    Vector3d ds;            // tangential increment
    double   Fsmag;         // magnitude of frictional/shear force
    Vector3d k;             // axis of rotation to rotate old normal to new normal
    double   sint, cost;
    int      ncontacts = 0;     // number of contact points
    Vector3d branchVec = _position - other._position;
    double GamaN = -2*sqrt(_kn*_mass*other._mass)/(_mass+other._mass)*log(_cresN)/sqrt(M_PI*M_PI+log(_cresN)*log(_cresN));
    penetration = branchVec.norm() - _radius - other._radius;
    auto search = _nodeFriction.find(other._id);
    if(penetration < 0.){
      /* normal force */
      normal = -branchVec / branchVec.norm();
      v = _velocity - other._velocity;
      Fn = penetration*normal*_kn - GamaN*normal.dot(v)*normal;
      /* tangential force */
      ds = (v - v.dot(normal)*normal) * dt;
      if(search != _nodeFriction.end()){
        Vector3d nodeNormal, nodeShear;
        std::tie(nodeNormal, nodeShear) = _nodeFriction[other._id];
        k = normal.cross(nodeNormal);
        sint = k.norm(); cost = sqrt(1-sint*sint); if(isnan(cost)){cost = 0.;}
        k = k/(sint+std::numeric_limits<double>::epsilon());
        newNodeShear = nodeShear*cost + k.cross(nodeShear)*sint + k*k.dot(nodeShear)*(1.0-cost);
      }
      newNodeNormal = normal;
      newNodeShear -= ds*_ks;
      newNodeContact = other._id;
      Fsmag = min(Fn.norm()*_mu, newNodeShear.norm());
      if(Fsmag > 0){
        Fs = Fsmag*newNodeShear/newNodeShear.norm(); newNodeShear = Fs;
      }
      force += (Fn + Fs);
      _nodeFriction[newNodeContact] = std::make_tuple(newNodeNormal, newNodeShear);
    }else if(search != _nodeFriction.end()){
      _nodeFriction.erase(other._id);
    }
    double ncmax = 16.;
    if((double)ncontacts > ncmax){
      force *= ncmax/(double)ncontacts;
    }
    _force += force;
    other._force -= force;
    stressVoigt(0) += force(0)*branchVec(0);
    stressVoigt(1) += force(1)*branchVec(1);
    stressVoigt(2) += force(2)*branchVec(2);
    stressVoigt(3) += 0.5*(force(1)*branchVec(2) + force(2)*branchVec(1));
    stressVoigt(4) += 0.5*(force(2)*branchVec(0) + force(0)*branchVec(2));
    stressVoigt(5) += 0.5*(force(1)*branchVec(0) + force(0)*branchVec(1));
  }

  void takeTimeStep(const double & gDamping, const double & dt){
    _velocity = 1./(1.+gDamping*dt/2.)*((1.-gDamping*dt/2.)*_velocity+dt*_force/_mass);
    _position += dt*_velocity;
    _force << 0.,0.,0.;
  }
  void applyAcceleration(const Vector3d & acceleration){
    _force += _mass*acceleration;
  }
  double computeKineticEnergy() const{
    double ke = 0.;
    ke += 0.5*_mass*_velocity.squaredNorm();
    ke = isnan(ke) ? 0. : ke;
    return ke;
  }
  void changeDensity(const double & density) {
  	_mass *= density/_density;
  	_density = density;
  }
  void changePosition(const Vector3d & position){
    _position = position;
  }
  void changeVelocity(const Vector3d & velocity){
    _velocity = velocity;
  }
  void changeId(const int & id){
  	_id = id;
  }
  void changeKn(const double & kn) {
  	_kn = kn;
  }
  void changeMu(const double & mu) {
  	_mu = mu;
  }
  void changeKs(const double & ks) {
  	_ks = ks;
  }
  void changeForce(const Vector3d & force){
    _force = force;
  }
  void changeMass(const double & mass){
    _mass = mass;
  }
  void changeFrictionTuple(const int & id, const std::tuple<Vector3d, Vector3d> & nodeFriction){
    _nodeFriction[id] = nodeFriction;
  }
  void eraseFrictionTuple(const int & id){
    _nodeFriction.erase(id);
  }
  const double & getMass() const {
  	return _mass;
  }
  const Vector3d & getPosition() const {
  	return _position;
  }
  const Vector3d & getVelocity() const {
  	return _velocity;
  }
  const double & getRadius() const {
  	return _radius;
  }
  const double & getKn() const {
  	return _kn;
  }
  const double & getKs() const {
  	return _ks;
  }
  const double & getMu() const {
  	return _mu;
  }
  const double & getDensity() const {
  	return _density;
  }
  const int & getId() const {
  	return _id;
  }
  const Vector3d & getForce() const{
    return _force;
  }
  const int & getNContact() const{
    return _nContact;
  }
  void getIfContact(int id, bool & found) const{
    auto search = _nodeFriction.find(id);
    found = ( search != _nodeFriction.end() );
  }
  std::tuple<Vector3d, Vector3d> getFrictionTuple(int id){
    return _nodeFriction[id];
  }
  void addForce(const Vector3d & df){
    _force += df;
  }

private:
  double 		       _mass;
	Vector3d 	       _position; 		// location of center of mass in real space
	Vector3d 	       _velocity;
	double 		       _radius;       // radius of bounding sphere
	double		       _kn; 		      // normal stiffness
	double		       _ks; 		      // shear stiffness (currently not used)
	double		       _mu; 		      // interparticle friction
	double 		       _density;      // density of particle (default = 1)
	int              _id;
  Vector3d         _force;        // force at the center of mass of a grain
  double           _cresN;
  int              _nContact;     //number of contact
  std::map<int, std::tuple<Vector3d, Vector3d> > _nodeFriction;
};
#endif /*GRAIN3D_H_*/
