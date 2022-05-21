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
#include "Levelset3d.h"

extern int openmpThreads;
extern int dynamicSchedule;


const bool debug = false;
class Grain3d{
public:
  Grain3d(){
    _radius   = 0;
    _rsq      = 0;
    _kn       = 0;
    _ks       = 0;
    _cresN    = 0;
    _mass     = 0;
    _mu       = 0;
    _id       = 0;
    _density  = 0;
    _lset     = Levelset3d();
    _quat     << 0,0,0,0;
    _omega    << 0,0,0;
    _momentInertia << 0,0,0;
    _position << 0,0,0;
    _velocity << 0,0,0;
    _force    << 0,0,0;
    _moment   << 0,0,0;
    _cmLset   << 0,0,0;
    _rotMatrix << 0.,0.,0.,0.,0.,0.,0.,0.,0.;
    _refPointList.clear();
    _refPointList.shrink_to_fit();
    _nodeShears.clear();
    _nodeShears.shrink_to_fit();
    _nodeContact.clear();
    _nodeContact.shrink_to_fit();
    _nodeNormals.clear();
    _nodeNormals.shrink_to_fit();
  }
  ~Grain3d(){
    _radius   = 0;
    _rsq      = 0;
    _kn       = 0;
    _ks       = 0;
    _cresN    = 0;
    _mass     = 0;
    _mu       = 0;
    _id       = 0;
    _density  = 0;
    _quat     << 0,0,0,0;
    _omega    << 0,0,0;
    _momentInertia << 0,0,0;
    _position << 0,0,0;
    _velocity << 0,0,0;
    _force    << 0,0,0;
    _moment   << 0,0,0;
    _cmLset   << 0,0,0;
    _rotMatrix << 0.,0.,0.,0.,0.,0.,0.,0.,0.;
    _refPointList.clear();
    _refPointList.shrink_to_fit();
    _nodeShears.clear();
    _nodeShears.shrink_to_fit();
    _nodeContact.clear();
    _nodeContact.shrink_to_fit();
    _nodeNormals.clear();
    _nodeNormals.shrink_to_fit();
  }
  void clearGrain(){
    _radius   = 0;
    _rsq      = 0;
    _kn       = 0;
    _ks       = 0;
    _cresN    = 0;
    _mass     = 0;
    _mu       = 0;
    _id       = 0;
    _density  = 0;
    _lset.clearLevelset();
    _quat     << 0,0,0,0;
    _omega    << 0,0,0;
    _momentInertia << 0,0,0;
    _position << 0,0,0;
    _velocity << 0,0,0;
    _force    << 0,0,0;
    _moment   << 0,0,0;
    _cmLset   << 0,0,0;
    _rotMatrix << 0.,0.,0.,0.,0.,0.,0.,0.,0.;
    _refPointList.clear();
    _refPointList.shrink_to_fit();
    _nodeShears.clear();
    _nodeShears.shrink_to_fit();
    _nodeContact.clear();
    _nodeContact.shrink_to_fit();
    _nodeNormals.clear();
    _nodeNormals.shrink_to_fit();
  }
  Grain3d(const Vector3d & position, const Vector4d & quat, const int& gid){
    _radius   = 0;
    _rsq      = 0;
    _kn       = 0;
    _ks       = 0;
    _cresN    = 0;
    _mass     = 0;
    _mu       = 0;
    _id       = gid;
    _density  = 0;
    _lset     = Levelset3d();
    _quat     = quat;
    _omega    << 0,0,0;
    _momentInertia << 0,0,0;
    _position = position;
    _velocity << 0,0,0;
    _force    << 0,0,0;
    _moment   << 0,0,0;
    _cmLset   << 0,0,0;
    _rotMatrix << 0.,0.,0.,
                  0.,0.,0.,
                  0.,0.,0.;
    _inertiaTensor << 0.,0.,0.,
                      0.,0.,0.,
                      0.,0.,0.;
    _refPointList.clear();
    _refPointList.shrink_to_fit();
    _nodeShears.clear();
    _nodeShears.shrink_to_fit();
    _nodeContact.clear();
    _nodeContact.shrink_to_fit();
    _nodeNormals.clear();
    _nodeNormals.shrink_to_fit();
    //_realGrain = false;
  }
  Grain3d(const double & mass, const Vector3d & position, const Vector3d & velocity,
          const Vector3d & momentInertia, const Vector4d & quat, const Vector3d & omega,
          const Vector3d & cmLset, const vector<Vector3d> & refPointList, const double & radius,
          const Levelset3d& lset, const double kn, const double ks, const double mu, const int & id):
          _mass(mass), _position(position), _velocity(velocity), _momentInertia(momentInertia),
          _quat(quat), _omega(omega), _cmLset(cmLset), _refPointList(refPointList), _radius(radius),
          _lset(lset), _kn(kn), _ks(ks), _mu(mu), _id(id){
    _cresN   = 0.6;
    _density = 1.0;
    _force   << 0.,0.,0.;
    _moment  << 0.,0.,0.;

    double maxR = 0;
    for(int i = 0; i < _refPointList.size(); ++i){
      if(_refPointList[i].norm()>maxR)
        maxR = _refPointList[i].norm();
    }
    _radius = maxR; _rsq = _radius * _radius;
    // get _rotMatrix from updateRotationMatrix(quat)
    _dOmega << 0, 0, 0;
    _rotMatrix << 0.,0.,0.,
                  0.,0.,0.,
                  0.,0.,0.;
    _inertiaTensor << _momentInertia(0), 0., 0., 0., _momentInertia(1), 0., 0., 0., _momentInertia(2);
    updateRotationMatrix(quat);
    // global angular velocity, and node coordinate
    _omegaGlobal = _rotMatrix * _omega;
    _pointList.resize(_refPointList.size());
    _nodeShears.resize(_pointList.size());
    _nodeContact.resize(_pointList.size());
    _nodeNormals.resize(_pointList.size());

    _nodeShears.assign(_nodeShears.size(),Vector3d(0,0,0));
    _nodeNormals.assign(_nodeNormals.size(),Vector3d(0,0,0));
    _nodeContact.assign(_nodeContact.size(),-1);

    for(int i = 0; i < _refPointList.size(); ++i){
      //_rotMatrix convert principal coordinate to current coordinate
      _pointList[i] = _rotMatrix*_refPointList[i]+_position;
    }
  }

  Grain3d(const double & mass, const Vector3d & position, const Vector3d & velocity, const double & radius, const double kn, const double mu, const int & id):
          _mass(mass), _position(position), _velocity(velocity), _radius(radius), _kn(kn), _mu(mu), _id(id){
    _cresN   = 0.6;
    _density = 1.0;
    _force   << 0.,0.,0.;
    _quat    << 0.,0.,0.,0.;
    _omega   << 0.,0.,0.;
    _rsq     = _radius*_radius;
  }

  Grain3d(const Grain3d& other){
    _mass = other._mass;
    _position = other._position;
    _velocity = other._velocity;
    _momentInertia = other._momentInertia;
    _quat = other._quat;
    _rotMatrix = other._rotMatrix;
    _inertiaTensor = other._inertiaTensor;
    _omega = other._omega;
    _dOmega = other._dOmega;
    _cmLset	= other._cmLset;
    _pointList = other._pointList;
    _refPointList = other._refPointList;
    _radius = other._radius;
    _rsq = other._rsq;
    _lset = other._lset;
    _kn = other._kn;
    _ks = other._ks;
    _mu = other._mu;
    _density = other._density;
    _omegaGlobal = other._omegaGlobal;
    _id = other._id;
    _force = other._force;
    _moment = other._moment;
    _nodeShears = other._nodeShears;
    _nodeContact = other._nodeContact;
    _nodeNormals = other._nodeNormals;
    _cresN = other._cresN;
    _pointList.shrink_to_fit();
    _refPointList.shrink_to_fit();
    _nodeShears.shrink_to_fit();
    _nodeContact.shrink_to_fit();
    _nodeNormals.shrink_to_fit();
  }
  //MAY NOT NEED IT, IT IS SLOW TO COPY
  Grain3d & operator=(const Grain3d& other){
    _mass = other._mass;
    _position = other._position;
    _velocity = other._velocity;
    _momentInertia = other._momentInertia;
    _quat = other._quat;
    _rotMatrix = other._rotMatrix;
    _inertiaTensor = other._inertiaTensor;
    _omega = other._omega;
    _dOmega = other._dOmega;
    _cmLset	= other._cmLset;
    _pointList = other._pointList; //vector<Vector3d>
    _refPointList = other._refPointList; //vector<Vector3d>
    _radius = other._radius;
    _rsq = other._rsq;
    _lset = other._lset;
    _kn = other._kn;
    _ks = other._ks;
    _mu = other._mu;
    _density = other._density;
    _omegaGlobal = other._omegaGlobal;
    _id = other._id;
    _force = other._force;
    _moment = other._moment;
    _nodeShears = other._nodeShears; //vector<Vector3d>
    _nodeContact = other._nodeContact; //vector<int>
    _nodeNormals = other._nodeNormals; //vector<Vector3d>
    _cresN = other._cresN;
    _pointList.shrink_to_fit();
    _refPointList.shrink_to_fit();
    _nodeShears.shrink_to_fit();
    _nodeContact.shrink_to_fit();
    _nodeNormals.shrink_to_fit();
    return *this;
  }
  //FUNCTIONS
  bool bCircleCheck(const Grain3d & other) const{
    Vector3d d = other.getPosition()-_position;
    double rsum = other.getRadius()+_radius;
    //if the center distance between two grains smaller than their sum of _radius
    //then two grains would be likely to contact, it is not guaranteed because _radius is overestimated
    //more check will follow
    return d.squaredNorm() < rsum*rsum;
  }

  // Checks contact between *this and other.
	// Compares points of *this to the level set of other.
	// If there is contact, updates force, which is the force on *this,
	// thisMoment, which is the moment on *this, and otherMoment, the moment on other.
  void findInterparticleForceMoment(Grain3d & other, const double & dt){
    //zero the ouputs force, thisMoment, otherMoment
    Vector3d force(0., 0., 0.);
    Vector3d thisMoment(0., 0., 0.);
    Vector3d otherMoment(0., 0., 0.);
    Vector3d ptThisCM;      // distance to the center of mass of *this in real space
    Vector3d ptOtherCM;     // distance to the center of mass of other in real space
    Vector3d ptOtherLset;   // point in the reference config of other's level set, should be in reeal/global frame
    double   penetration;   // penetration amount
    Vector3d normal;        // surface normal pointing out of other in the reference config, point toward *this
    Vector3d Fn;            // normal force form a single point
    Vector3d v;             // velocity of this grain wrt to other grain
    Vector3d Fs;            // vector of frictional/shear force
    Vector3d ds;            // tangential increment
    double   Fsmag;         // magnitude of frictional/shear force
    Vector3d k;             // axis of rotation to rotate old normal to new normal
    double   sint, cost;
    int      ncontacts = 0; // number of contact points

    double GamaN = -2*sqrt(_kn*_mass*other.getMass()/(_mass+other.getMass()))*log(_cresN)/sqrt(M_PI*M_PI+log(_cresN)*log(_cresN));
    for(int i = 0; i < _pointList.size(); i+=5){
      // all _position and _pointList are in real/global frame
      ptOtherCM = _pointList[i]-other.getPosition();
      if(ptOtherCM.squaredNorm()<other.getRsq()){
        //rotate back to principal frame of other
        ptOtherLset = other.getRotMatrix().transpose()*ptOtherCM + other.getCmLset();
        //find amount of penetration, normal: point out of other, point to *this
        if(other.getLset().findPenetration(ptOtherLset, penetration, normal)){
          //return yes, then ptThisCM distance to the center of *this in real frame
          ptThisCM = _pointList[i]-_position; //checkflag = true;
          //rotate direction into real frame
          normal = -other.getRotMatrix()*normal;
          //relative velocity between two grains
          v = _velocity-other.getVelocity()+_omegaGlobal.cross(ptThisCM)-other.getOmegaGlobal().cross(ptOtherCM);
          Fn = penetration*normal*_kn - GamaN*normal.dot(v)*normal;
          //cout<<" gid: "<<_id<<" incremental normal force is: "<<Fn(0)<<" "<<Fn(1)<<" "<<Fn(2)<<endl;
          //Fn in real frame, force, thisMoment, otherMoment in real frame
          force += Fn;
          thisMoment  += ptThisCM.cross(Fn);
          otherMoment += ptOtherCM.cross(-Fn);
          ncontacts ++;
          //use dt to find the tangential incremental displacement
          ds = (v - v.dot(normal)*normal)*dt;
          //cout<<" gid: "<<_id<<" incremental shear displacement is: "<<ds(0)<<" "<<ds(1)<<" "<<ds(2)<<endl;
          if(_nodeNormals[i].norm()>0){
            k = normal.cross(_nodeNormals[i]);
            sint = k.norm();
            cost = sqrt(1-sint*sint);
            if(isnan(cost)){cost = 0.;}
            k = k/(sint+std::numeric_limits<double>::epsilon());
            // this is Rodrigues' rotation formula
            _nodeShears[i] = _nodeShears[i]*cost+k.cross(_nodeShears[i])*sint+k*k.dot(_nodeShears[i])*(1.0-cost);
          }
          _nodeNormals[i] = normal;        // normal direction of the ith node in real frame
          _nodeShears[i] -= ds*_ks;        // shear force of the ith node in real frame
          _nodeContact[i] = other.getId(); //contactId of the ith node
          //the magnitude of shear force cannot be larger than mu*Fn
          Fsmag = min(Fn.norm()*_mu, _nodeShears[i].norm());
          if(Fsmag > 0){
            Fs = Fsmag*_nodeShears[i]/_nodeShears[i].norm();
            //update _nodeShears, update friction force to force, thisMoment, otherMoment
            _nodeShears[i] = Fs;
            force += Fs;
            thisMoment += ptThisCM.cross(Fs);
            otherMoment += ptOtherCM.cross(-Fs);
          }
        }
      }
    }
    // prevent instability by limiting effective number of contact points based on mass
    double ncmax = 16.;
    if((double)ncontacts > ncmax){
      force *= ncmax/(double)ncontacts;
      thisMoment *= ncmax/(double)ncontacts;
      otherMoment *= ncmax/(double)ncontacts;
    }
    force *= 5;
    thisMoment *= 5;
    otherMoment *= 5;
    _force += force;
    other._force -= force;
    _moment += thisMoment;
    other._moment += otherMoment;
    if(isnan(abs(_force(0))) || isnan(abs(_force(1))) || isnan(abs(_force(2)))){
      cout<<"gid: "<<_id<<" has NAN force."<<endl;
      _force << 0., 0., 0.;
      _moment << 0., 0., 0.;
      other._force << 0., 0., 0.;
      other._moment << 0., 0., 0.;
    }
  }

  void takeTimestep(const double & gDamping, const double & dt){
    _velocity = 1./(1.+gDamping*dt/2.)*((1.-gDamping*dt/2.)*_velocity+dt*_force/_mass);
    Vector3d principleMoment = _rotMatrix.transpose()*_moment;
    Vector3d omegaN = _omega+_dOmega/2;
    for (int i = 0; i < 5; i++) {
  		_dOmega(0) = (principleMoment(0) + omegaN(1)*omegaN(2)*(_momentInertia(1)-_momentInertia(2)) - gDamping*_momentInertia(0)*omegaN(0) )*dt/_momentInertia(0);
  		_dOmega(1) = (principleMoment(1) + omegaN(2)*omegaN(0)*(_momentInertia(2)-_momentInertia(0)) - gDamping*_momentInertia(1)*omegaN(1) )*dt/_momentInertia(1);
  		_dOmega(2) = (principleMoment(2) + omegaN(0)*omegaN(1)*(_momentInertia(0)-_momentInertia(1)) - gDamping*_momentInertia(2)*omegaN(2) )*dt/_momentInertia(2);
  		omegaN = _omega + _dOmega/2.;
    }

    if(isinf(_dOmega(0)) || isnan(_dOmega(0)) || isinf(_dOmega(1)) || isnan(_dOmega(1)) || isinf(_dOmega(2)) || isnan(_dOmega(2))){
      if(debug){
        cout << "oops... bad id: " << _id << " goes into infinity..." << endl;
        cout << "omega is: "<<_omega(0)<<" "<<_omega(1)<<" "<<_omega(2)<<endl;
        cout << "dOmega is: "<<_dOmega(0)<<" "<<_dOmega(1)<<" "<<_dOmega(2)<<endl;
        cout << "velocity is: "<<_velocity(0)<<" "<<_velocity(1)<<" "<<_velocity(2)<<endl;
        cout << "gDamping is: "<<gDamping<<" force is: "<<_force(0)<<" "<<_force(1)<<" "<<_force(2)<<" mass is: "<<_mass<<endl;
        cout<<"############## END OUTPUT ################"<<endl;
      }
      _dOmega << 0.,0.,0.;
      _velocity << 0.,0.,0.;
      _omega << 0.,0.,0.;
      _quat << 0.,0.,0.,1.;
    }
    _omega += _dOmega;
  	// everything after this should be correct
  	// update quaternions and rotation matrix
  	// convert quat to go from global to reference (the opposite of what it usually is, but that's needed for this bit here)
  	_quat << -_quat(0), -_quat(1), -_quat(2), _quat(3) ;
  	double Bx = dt/4.*_omega(0);
  	double By = dt/4.*_omega(1);
  	double Bz = dt/4.*_omega(2);
  	double c1 =  _quat(0) + Bz*_quat(1) - Bx*_quat(2) - By*_quat(3);
  	double c2 = -Bz*_quat(0) + _quat(1) - By*_quat(2) + Bx*_quat(3);
  	double c3 =  Bx*_quat(0) + By*_quat(1) + _quat(2) + Bz*_quat(3);
  	double c4 =  By*_quat(0) - Bx*_quat(1) - Bz*_quat(2) + _quat(3);
  	double detB = 1 + 2*Bx*Bx + 2*By*By + 2*Bz*Bz + 2*Bx*Bx*By*By + 2*By*By*Bz*Bz + 2*Bx*Bx*Bz*Bz + Bx*Bx*Bx*Bx + By*By*By*By + Bz*Bz*Bz*Bz;
  	double bfac = (1 + Bx*Bx + By*By + Bz*Bz)/detB;
  	_quat(0) = ( c1 + c2*Bz - c3*Bx - c4*By)*bfac;
  	_quat(1) = (-c1*Bz + c2 - c3*By + c4*Bx)*bfac;
  	_quat(2) = ( c1*Bx + c2*By + c3 + c4*Bz)*bfac;
  	_quat(3) = ( c1*By - c2*Bx - c3*Bz + c4)*bfac;
  	// normalize _quat and use it to update rotation matrix
  	_quat = _quat/_quat.norm();
  	// convert quaternion back (going from reference to global, which is what it usually is)
  	_quat << -_quat(0), -_quat(1), -_quat(2), _quat(3);
    updateRotationMatrix(_quat);
    _omegaGlobal = _rotMatrix*_omega;
    _position += dt*_velocity;
  }

  void applyAcceleration(const Vector3d & acceleration){
    _force += _mass*acceleration;
  }
  // this is because _position and _rotMatrix changed

  void updatePoints(){
    for(int i = 0; i < _pointList.size(); ++i){
      _pointList[i] = _rotMatrix*_refPointList[i]+_position;
    }
  }
  //this is only one needed in simple simulation

  double computeKineticEnergy() const{
    double ke = 0;
    ke += 0.5*_mass*_velocity.squaredNorm();
    for(int i = 0; i < 3; ++i){
      ke += 0.5*_omega(i)*_momentInertia(i)*_omega(i);
    }
    return ke;
  }
  void changeDensity(const double & density) {
  	_mass *= density/_density;
  	_momentInertia *= density/_density;
  	_density = density;
  }
  void changeRotation(const Vector4d & quat) {
  	_quat = quat/quat.norm();
  	updateRotationMatrix(_quat);
  }
  void changeRotation(const double quat[]) {
    _quat << quat[0],quat[1],quat[2],quat[3];
  	_quat = _quat/_quat.norm();
  	updateRotationMatrix(_quat);
  }
  void changePosition(const Vector3d & position){
    _position = position;
  }
  void changeQuat(const Vector4d & quat){
    _quat = quat;
  }
  void changeQuat(const double quat[]){ //this is used because MPI communication cannot pass Vector4d
    for(size_t i = 0; i < 4; ++i) _quat(i) = quat[i];
  }
  void changeVelocity(const Vector3d & velocity){
    _velocity = velocity;
  }
  void changeOmega(const Vector3d & omega){
    _omega = omega;
    _omegaGlobal = _rotMatrix * _omega;
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
  void changeMoment(const Vector3d & moment){
    _moment = moment;
  }
  void changeDOmega(const Vector3d & dOmega){
    _dOmega = dOmega;
  }
  void addForce(const Vector3d & force){
    _force += force;
  }
  void addMoment(const Vector3d & moment){
    _moment += moment;
  }
  void changeShearHist(const vector<Vector3d> & nodeShears, const vector<int> & nodeContact, const vector<Vector3d> & nodeNormals){
    _nodeShears = nodeShears;
    _nodeContact = nodeContact;
    _nodeNormals = nodeNormals;
  }

  void clearFriction() {
    _nodeShears.assign(_nodeShears.size(),Vector3d(0,0,0));
    _nodeNormals.assign(_nodeNormals.size(),Vector3d(0,0,0));
    _nodeContact.assign(_nodeContact.size(),-1);
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
  const Vector4d & getQuat() const {
  	return _quat;
  }
  const Matrix3d & getRotMatrix() const {
  	return _rotMatrix;
  }
  const Vector3d & getOmega() const {
  	return _omega;
  }
  const Vector3d & getOmegaGlobal() const {
  	return _omegaGlobal;
  }
  const Vector3d & getMomentInertia() const {
  	return _momentInertia;
  }
  const Vector3d & getCmLset() const {
  	return _cmLset;
  }
  const double & getRadius() const {
  	return _radius;
  }
  const double & getRsq() const {
  	return _rsq;
  }
  const Levelset3d getLset() const {
  	return _lset;
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
  const vector<Vector3d> & getPointList() const {
  	return _pointList;
  }
  const int & getId() const {
  	return _id;
  }
  const vector<int> & getNodeContact(){
  	return _nodeContact;
  }
  vector<int> & getNodeContactNonConst(){
  	return _nodeContact;
  }
  const vector<Vector3d> & getNodeShears(){
  	return _nodeShears;
  }
  vector<Vector3d> & getNodeShearsNonConst(){
  	return _nodeShears;
  }
  const vector<Vector3d> & getNodeNormals(){
  	return _nodeNormals;
  }
  vector<Vector3d> & getNodeNormalsNonConst(){
  	return _nodeNormals;
  }
  const Vector3d & getForce() const{
    return _force;
  }
  const Vector3d & getMoment() const{
    return _moment;
  }
  void setRealGrain(bool real){
    _realGrain = real;
  }
  bool getRealGrain(){
    return _realGrain;
  }
//may be useful for eigen alignment issue. but -std=c++17 would solve everything
EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
  void updateQuat(const Matrix3d & R) {
    double tr = R.trace();
    _quat(3) = sqrt(tr+1.)/2.;
    _quat(0) = (R(0,2)-R(2,0))/(4*_quat(3));
    _quat(1) = (R(1,2)-R(2,1))/(4*_quat(3));
    _quat(2) = (R(0,1)-R(1,0))/(4*_quat(3));
  }
  void updateRotationMatrix(const Vector4d & quat) {
    _rotMatrix(0,0) = -quat(0)*quat(0) + quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
    _rotMatrix(0,1) = -2*(quat(0)*quat(1) - quat(2)*quat(3));
    _rotMatrix(0,2) =  2*(quat(1)*quat(2) + quat(0)*quat(3));
    _rotMatrix(1,0) = -2*(quat(0)*quat(1) + quat(2)*quat(3));
    _rotMatrix(1,1) =  quat(0)*quat(0) - quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
    _rotMatrix(1,2) = -2*(quat(0)*quat(2) - quat(1)*quat(3));
    _rotMatrix(2,0) =  2*(quat(1)*quat(2) - quat(0)*quat(3));
    _rotMatrix(2,1) = -2*(quat(0)*quat(2) + quat(1)*quat(3));
    _rotMatrix(2,2) = -quat(0)*quat(0) - quat(1)*quat(1) + quat(2)*quat(2) + quat(3)*quat(3);
  }
  Vector4d 	       _quat; 			    // quaternion representing the rotation in (X,Y,Z,W) form; principle moment of inertia and the rotation from principal -> current respectively
  double 		       _mass;
	Vector3d 	       _position; 		  // location of center of mass in real space
	Vector3d 	       _velocity;
	Vector3d 	       _momentInertia;  // moment of inertia in principal frame (purely diagonal terms), used for energy computation
  Matrix3d         _inertiaTensor;  // inertial tensor in global
	Matrix3d		     _rotMatrix; 	    // rotation matrix from principle frame to current frame
	Vector3d		     _omega; 			    // angular velocity IN THE PRINCIPLE FRAME
	Vector3d		     _dOmega; 		    // change in angular velocity between timesteps (needed for predictor-corrector algorithm) IN THE PRINCIPLE FRAME
	Vector3d		     _cmLset; 		    // center of mass wrt the level set reference configuration
	vector<Vector3d> _pointList; 		  // list of points comprising the grain in real space (translated and rotated)
	vector<Vector3d> _refPointList;   // list of points in reference config (center of mass is at (0,0,0) and I is diagonal)
	double 		       _radius;         // radius of bounding sphere
	double		       _rsq;		        // squared radius
	Levelset3d       _lset;	          // level set of grain
	double		       _kn; 		        // normal stiffness
	double		       _ks; 		        // shear stiffness (currently not used)
	double		       _mu; 		        // interparticle friction
	double 		       _density;        // density of particle (default = 1)
	Vector3d         _omegaGlobal;
	int              _id;
  Vector3d         _force;          // force at the center of mass of a grain
  Vector3d         _moment;         // moment at the center of mass of a grain
	vector<Vector3d> _nodeShears;     // shear force at each node
	vector<int>      _nodeContact;    // index of grain the node is contacting
	vector<Vector3d> _nodeNormals;    // contact normals at each node
  double           _cresN;
  bool             _realGrain;
};
#endif /*GRAIN3D_H_*/
