#ifndef GRAIN3D_HPP_
#define GRAIN3D_HPP_
#include "Levelset3d.hpp"
#include "utilities.hpp"

struct ContactInfo;

class Grain3d{
public:
  Grain3d();
  Grain3d(const Vector3d & position, const Vector4d & quat, const int& gid);
  Grain3d(const double & mass, const double & density, const Vector3d & position, const Vector3d & velocity,
          const Vector3d & momentInertia, const Vector4d & quat, const Vector3d & omega, const double & mu, const int & id);
  Grain3d(const double & mass, const double & density, const Vector3d & position, const Vector3d & velocity,
          const Vector3d & momentInertia, const Vector4d & quat, const Vector3d & omega,
          const Vector3d & cmLset, const vector<Vector3d> & refPointList, const double & radius,
          const Levelset3d& lset, const double mu, const int & id);
  Grain3d(const Grain3d& other);
  Grain3d & operator=(const Grain3d& other);
  bool bCircleCheck(const Grain3d & other) const;

  void findInterGrainContact(const Grain3d & other, vector<ContactInfo> & contactList, vector<int> & contactId);
  void takeTimestep(const double & dt);
  void takeTempTimestep(const double & dt);
  void updatePoints();
  void updateTempPoints();
  double computeKineticEnergy();
  void changeMass(const double & mass);
  void changePosition(const Vector3d & position);
  void changeVelocity(const Vector3d & velocity);
  void changeOmega(const Vector3d & omega);
  void changeTempOmega(const Vector3d & omega);
  void changeOmegaGlobal(const Vector3d & omegaGlobal);
  void changeDensity(const double & density);
  void changeRotation(const Vector4d & quat);
  void changeRotation(const double quat[]);
  void changeMu(const double & mu);
  void changeQuat(const Vector4d & quat);
  void changeQuat(const double quat[]);
  void addVelocityCompression(const Vector3d & vel, const Vector3d & normal, const int & size);
  void addVelocitySeparation(const Vector3d & vel);
  void addVelocity(const Vector3d & vel);
  void addOmegaGlobal(const Vector3d & omg);
  void addVelocityAcceleration(const Vector3d & acc, double dt);
  const double   & getRadius() const;
  const Vector4d & getQuat() const;
  const Vector3d & getPosition() const;
  const Vector3d & getTempPosition() const;
  const Vector3d & getVelocity() const;
  const Vector3d & getOmega() const;
  const double   & getRsq() const;
  const Vector3d & getCmLset() const;
  const Matrix3d & getRotMatrix() const;
  const Matrix3d & getTempRotMatrix() const;
  const Matrix3d & getInertiaTensor() const;
  const Vector3d & getMomentInertia() const;
  const Matrix3d   getGlobalI() const;
  const double   & getMass() const;
  const int      & getId() const;
  const Vector3d & getOmegaGlobal() const;
  const Levelset3d & getLset() const;
  const vector<Vector3d> & getPointList() const;
  void setRealGrain(bool realGrain);
  const bool & isRealGrain() const;

private:
  void updateQuat(const Matrix3d & R);
  void updateRotationMatrix(const Vector4d & quat);
  void updateTempRotationMatrix(const Vector4d & quat);
  int              _id;
  Vector4d 	       _quat; 			  // quaternion representing the rotation in (X,Y,Z,W) form; principle moment of inertia and the rotation from principal -> current respectively
  double 		       _mass;
	Vector3d 	       _position; 		// location of center of mass in real space
	Vector3d 	       _velocity;
	Vector3d 	       _momentInertia;// moment of inertia in principal frame (purely diagonal terms), used for energy computation
  Matrix3d         _inertiaTensor;// inertial tensor in global
	Matrix3d		     _rotMatrix; 	  // rotation matrix from principle frame to current frame
	Vector3d		     _omega; 			  // angular velocity IN THE PRINCIPLE FRAME
	Vector3d		     _cmLset; 		  // center of mass wrt the level set reference configuration
	vector<Vector3d> _pointList; 		// list of points comprising the grain in real space (translated and rotated)
	vector<Vector3d> _refPointList; // list of points in reference config (center of mass is at (0,0,0) and I is diagonal)
	double 		       _radius;       // radius of bounding sphere
	double		       _rsq;		      // squared radius
	Levelset3d       _lset;	        // level set of grain
	double		       _mu; 		      // interparticle friction
	double 		       _density;      // density of particle (default = 1)
	Vector3d         _omegaGlobal;
  double           _cresN;
  Vector4d         _tempQuat;
  Vector3d         _tempPosition;
  Matrix3d         _tempRotMatrix;
  bool             _realGrain;
};

Grain3d::Grain3d(){
  _radius   = 0;
  _rsq      = 0;
  _mass     = 0;
  _mu       = 0;
  _id       = 0;
  _density  = 1.;
  _lset     = Levelset3d();
  _quat     << 0,0,0,0;
  _omega    << 0,0,0;
  _momentInertia << 0,0,0;
  _position << 0,0,0;
  _velocity << 0,0,0;
  _cmLset   << 0,0,0;
  _rotMatrix << 1.,0.,0.,
                0.,1.,0.,
                0.,0.,1.;
  _inertiaTensor << 0.,0.,0.,
                    0.,0.,0.,
                    0.,0.,0.;
  _refPointList.clear();
  _refPointList.shrink_to_fit();
}

Grain3d::Grain3d(const Vector3d & position, const Vector4d & quat, const int& gid){
  _radius   = 0;
  _rsq      = 0;
  _mass     = 0;
  _mu       = 0;
  _id       = gid;
  _density  = 1.;
  _lset     = Levelset3d();
  _quat     = quat;
  _omega    << 0,0,0;
  _momentInertia << 0,0,0;
  _position = position;
  _velocity << 0,0,0;
  _cmLset   << 0,0,0;
  _rotMatrix << 1.,0.,0.,
                0.,1.,0.,
                0.,0.,1.;
  _inertiaTensor << 0.,0.,0.,
                    0.,0.,0.,
                    0.,0.,0.;
  _refPointList.clear();
  _refPointList.shrink_to_fit();
}

Grain3d::Grain3d(const double & mass, const double & density, const Vector3d & position, const Vector3d & velocity,
        const Vector3d & momentInertia, const Vector4d & quat, const Vector3d & omega, const double & mu, const int & id):
        _mass(mass), _position(position), _velocity(velocity), _momentInertia(momentInertia),
        _quat(quat), _omega(omega), _mu(mu), _id(id){
          _rotMatrix << 1.,0.,0.,
                        0.,1.,0.,
                        0.,0.,1.;
          _density = 1.;
          changeDensity(density);
          _inertiaTensor << _momentInertia(0), 0., 0., 0., _momentInertia(1), 0., 0., 0., _momentInertia(2);
          updateRotationMatrix(quat);
          _omegaGlobal = _rotMatrix * _omega;
        }

Grain3d::Grain3d(const double & mass, double const & density, const Vector3d & position, const Vector3d & velocity,
        const Vector3d & momentInertia, const Vector4d & quat, const Vector3d & omega,
        const Vector3d & cmLset, const vector<Vector3d> & refPointList, const double & radius,
        const Levelset3d& lset, const double mu, const int & id):
        _mass(mass), _position(position), _velocity(velocity), _momentInertia(momentInertia),
        _quat(quat), _omega(omega), _cmLset(cmLset), _refPointList(refPointList), _radius(radius),
        _lset(lset), _mu(mu), _id(id){
  double maxR = 0;
  for(int i = 0; i < _refPointList.size(); ++i){
    if(_refPointList[i].norm()>maxR)
      maxR = _refPointList[i].norm();
  }
  _radius = maxR;
  _rsq = _radius * _radius;
  _rotMatrix << 1.,0.,0.,
                0.,1.,0.,
                0.,0.,1.;
  _density = 1.;
  changeDensity(density);
  _inertiaTensor << _momentInertia(0), 0., 0., 0., _momentInertia(1), 0., 0., 0., _momentInertia(2);
  updateRotationMatrix(quat);
  _tempPosition = _position;
  _tempQuat = _quat;
  _tempRotMatrix = _rotMatrix;
  _omegaGlobal = _rotMatrix * _omega;
  _pointList.resize(_refPointList.size());

  for(int i = 0; i < _refPointList.size(); ++i){
    _pointList[i] = _rotMatrix*_refPointList[i]+_position;
  }
}

Grain3d::Grain3d(const Grain3d& other){
  _mass = other._mass;
  _position = other._position;
  _velocity = other._velocity;
  _momentInertia = other._momentInertia;
  _quat = other._quat;
  _rotMatrix = other._rotMatrix;
  _inertiaTensor = other._inertiaTensor;
  _omega = other._omega;
  _cmLset	= other._cmLset;
  _pointList = other._pointList;
  _refPointList = other._refPointList;
  _radius = other._radius;
  _rsq = other._rsq;
  _lset = other._lset;
  _mu = other._mu;
  _density = other._density;
  _omegaGlobal = other._omegaGlobal;
  _id = other._id;
  _pointList.shrink_to_fit();
  _refPointList.shrink_to_fit();
}

Grain3d & Grain3d::operator=(const Grain3d& other){
  _mass = other._mass;
  _position = other._position;
  _velocity = other._velocity;
  _momentInertia = other._momentInertia;
  _quat = other._quat;
  _rotMatrix = other._rotMatrix;
  _inertiaTensor = other._inertiaTensor;
  _omega = other._omega;
  _cmLset	= other._cmLset;
  _pointList = other._pointList;
  _refPointList = other._refPointList;
  _radius = other._radius;
  _rsq = other._rsq;
  _lset = other._lset;
  _mu = other._mu;
  _density = other._density;
  _omegaGlobal = other._omegaGlobal;
  _id = other._id;
  _pointList.shrink_to_fit();
  _refPointList.shrink_to_fit();
  return *this;
}

bool Grain3d::bCircleCheck(const Grain3d & other) const{
  Vector3d d = other.getPosition()-_position;
  double rsum = other.getRadius()+_radius;
  return d.squaredNorm() < rsum*rsum;
}

void Grain3d::findInterGrainContact(const Grain3d & other, vector<ContactInfo> & contactList, vector<int> & contactId ){
  Vector3d ptThisCM;      // distance to the center of mass of *this in real space
  Vector3d ptOtherCM;     // distance to the center of mass of other in real space
  Vector3d ptOtherLset;   // point in the reference config of other's level set, should be in reeal/global frame
  double   penetration = 0.;   // penetration amount
  Vector3d normal;        // surface normal pointing out of other in the reference config, point toward *this
  Vector3d velocity;      // velocity of this grain wrt to other grain
  vector<int> contactNumbers;
  vector<Vector3d> contactNormals;
  vector<double> contactPenetration;
  vector<double> elasticEnergy;
  vector<Vector3d> branchVector1; //ptThisCM;
  vector<Vector3d> branchVector2; //ptOtherCM;
  for(int i = 0; i < _pointList.size(); ++i){
    ptOtherCM = _pointList[i]-other.getPosition();
    if(ptOtherCM.squaredNorm()<other.getRsq()){
      ptOtherLset = other.getRotMatrix().transpose()*ptOtherCM + other.getCmLset();
      if( other.getLset().findPenetration(ptOtherLset, penetration, normal) ){
        ptThisCM = _pointList[i]-_position;
        normal   = other.getRotMatrix()*normal;
        /*case 1, no contact yet*/
        if(contactNumbers.size() == 0){
          contactNumbers.push_back(1);
          contactNormals.push_back(normal);
          contactPenetration.push_back(penetration);
          elasticEnergy.push_back(0.);
          branchVector1.push_back(ptThisCM);
          branchVector2.push_back(ptOtherCM);
        }else{
          /*case 2, ptThisCM is similar*/
          bool flag = false;
          for(int j = 0; j < contactNumbers.size(); ++j){
            Vector3d r = branchVector1[j]/contactNumbers[j];
            double similarity = r.dot(ptThisCM)/r.norm()/ptThisCM.norm();
            if(similarity > 0.9){
              contactNumbers[j] += 1;
              contactNormals[j] += normal;
              contactPenetration[j] += penetration;
              branchVector1[j]  += ptThisCM;
              branchVector2[j]  += ptOtherCM;
              flag = true;
              break;
            }
          }
          /*case 3, new contact*/
          if(!flag){
            contactNumbers.push_back(1);
            contactNormals.push_back(normal);
            contactPenetration.push_back(penetration);
            elasticEnergy.push_back(0.);
            branchVector1.push_back(ptThisCM);
            branchVector2.push_back(ptOtherCM);
          }
        }
      }
    }
    penetration = 0.;
  }
  /* check convexity */
  for(int i = 0; i < contactNumbers.size(); ++i){
    contactNormals[i] = contactNormals[i]/contactNormals[i].norm();
  }
  for(int i = 0; i < contactNumbers.size(); ++i){
    for(int j = i+1; j < contactNumbers.size(); ++j){
      double similarity = contactNormals[i].dot(contactNormals[j]);
      if(similarity < similarityThreshold){
        return;
      }
    }
  }
  bool idFlag = false;
  for(int i = 0; i < contactNumbers.size(); i+=ceil( contactNumbers.size()/5. ) ){
    Vector3d bv1 = branchVector1[i]/contactNumbers[i];
    Vector3d bv2 = branchVector2[i]/contactNumbers[i];
    Vector3d n   = contactNormals[i];
    double avgPenetration = contactPenetration[i]/contactNumbers[i];
    velocity = this->getVelocity()-other.getVelocity()+this->getOmegaGlobal().cross(bv1)-other.getOmegaGlobal().cross(bv2);
    if(n.dot(velocity) < -threshold){
      idFlag = true;
      if(_mass < other._mass){
        contactList.emplace_back(_id, other._id, bv1, bv2, n, velocity, n.dot(velocity), elasticEnergy[i]);
      }else{
        contactList.emplace_back(other._id, _id, bv2, bv1, -n, -velocity, n.dot(velocity), elasticEnergy[i]);
      }
    }
  }
  if(idFlag){
    contactId.push_back(other._id);
  }
}

void Grain3d::takeTimestep(const double & dt){
  if(isnan(_velocity.norm()) || isnan(_omega.norm())){
    cout<<"Grain "<<_id<<" goes into infinity."<<endl;
    _velocity << 0.,0.,0.;
    _omega << 0.,0.,0.;
    _omegaGlobal << 0.,0.,0.;
  }
  _omega = _rotMatrix.transpose()*_omegaGlobal;
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
	_quat = _quat/_quat.norm();
	_quat << -_quat(0), -_quat(1), -_quat(2), _quat(3);
  _position += dt*_velocity;
  updateRotationMatrix(_quat);
  updatePoints(); // may not need it
}

void Grain3d::takeTempTimestep(const double & dt){
  if(isnan(_velocity.norm()) || isnan(_omega.norm())){
    cout<<"Grain "<<_id<<" goes into infinity."<<endl;
    _velocity << 0.,0.,0.;
    _omega << 0.,0.,0.;
    _omegaGlobal << 0.,0.,0.;
  }
  _tempQuat = _quat; _tempPosition = _position;
  _omega = _rotMatrix.transpose()*_omegaGlobal;
	_tempQuat << -_tempQuat(0), -_tempQuat(1), -_tempQuat(2), _tempQuat(3) ;
	double Bx = dt/4.*_omega(0);
	double By = dt/4.*_omega(1);
	double Bz = dt/4.*_omega(2);
	double c1 =  _tempQuat(0) + Bz*_tempQuat(1) - Bx*_tempQuat(2) - By*_tempQuat(3);
	double c2 = -Bz*_tempQuat(0) + _tempQuat(1) - By*_tempQuat(2) + Bx*_tempQuat(3);
	double c3 =  Bx*_tempQuat(0) + By*_tempQuat(1) + _tempQuat(2) + Bz*_tempQuat(3);
	double c4 =  By*_tempQuat(0) - Bx*_tempQuat(1) - Bz*_tempQuat(2) + _tempQuat(3);
	double detB = 1 + 2*Bx*Bx + 2*By*By + 2*Bz*Bz + 2*Bx*Bx*By*By + 2*By*By*Bz*Bz + 2*Bx*Bx*Bz*Bz + Bx*Bx*Bx*Bx + By*By*By*By + Bz*Bz*Bz*Bz;
	double bfac = (1 + Bx*Bx + By*By + Bz*Bz)/detB;
	_tempQuat(0) = ( c1 + c2*Bz - c3*Bx - c4*By)*bfac;
	_tempQuat(1) = (-c1*Bz + c2 - c3*By + c4*Bx)*bfac;
	_tempQuat(2) = ( c1*Bx + c2*By + c3 + c4*Bz)*bfac;
	_tempQuat(3) = ( c1*By - c2*Bx - c3*Bz + c4)*bfac;
	_tempQuat = _tempQuat/_tempQuat.norm();
	_tempQuat << -_tempQuat(0), -_tempQuat(1), -_tempQuat(2), _tempQuat(3);
  _tempPosition += dt*_velocity;
  updateTempRotationMatrix(_tempQuat);
  updateTempPoints();
}

void Grain3d::updatePoints(){
  for(int i = 0; i < _pointList.size(); ++i){
    _pointList[i] = _rotMatrix*_refPointList[i]+_position;
  }
}

void Grain3d::updateTempPoints(){
  for(int i = 0; i < _pointList.size(); ++i){
    _pointList[i] = _tempRotMatrix*_refPointList[i]+_tempPosition;
  }
}

double Grain3d::computeKineticEnergy(){
  double ke = 0;
  ke += 0.5*_mass*_velocity.squaredNorm();
  for(int i = 0; i < 3; ++i){
    ke += 0.5*_omega(i)*_momentInertia(i)*_omega(i);
  }
  if(isnan(ke) || isnan(-ke)){
    //cout<<"Grain "<<_id<<" goes into infinity."<<endl;
    ke = 0;
    _velocity << 0.,0.,0.;
    _omega << 0.,0.,0.;
    _omegaGlobal << 0.,0.,0.;
  }
  return ke;
}

void Grain3d::changeMass(const double & mass){
  _mass = mass;
}

void Grain3d::changePosition(const Vector3d & position){
  _position = _position;
}

void Grain3d::changeVelocity(const Vector3d & velocity){
  _velocity = velocity;
}

void Grain3d::changeOmega(const Vector3d & omega){
  _omega = omega;
  _omegaGlobal = _rotMatrix * _omega;
}

void Grain3d::changeTempOmega(const Vector3d & omega){
  _omega = omega;
  _omegaGlobal = _tempRotMatrix * _omega;
}

void Grain3d::changeOmegaGlobal(const Vector3d & omegaGlobal){
  _omegaGlobal = omegaGlobal;
  _omega = _rotMatrix.transpose() * _omegaGlobal;
}

void Grain3d::changeDensity(const double & density) {
	_mass *= density/_density;
	_momentInertia *= density/_density;
	_density = density;
}

void Grain3d::changeRotation(const Vector4d & quat) {
	_quat = quat/quat.norm();
	updateRotationMatrix(_quat);
  _omegaGlobal = _rotMatrix * _omega;
}

void Grain3d::changeRotation(const double quat[]) {
  _quat << quat[0],quat[1],quat[2],quat[3];
	_quat = _quat/_quat.norm();
	updateRotationMatrix(_quat);
  _omegaGlobal = _rotMatrix * _omega;
}

const double   & Grain3d::getRadius() const{
  return _radius;
}

const Vector4d & Grain3d::getQuat() const{
  return _quat;
}

const Vector3d & Grain3d::getPosition() const{
  return _position;
}

const Vector3d & Grain3d::getTempPosition() const{
  return _tempPosition;
}

const Vector3d & Grain3d::getVelocity() const{
  return _velocity;
}

const Vector3d & Grain3d::getOmega() const{
  return _omega;
}

const double   & Grain3d::getRsq() const{
  return _rsq;
}

const Vector3d & Grain3d::getCmLset() const{
  return _cmLset;
}

const Matrix3d & Grain3d::getRotMatrix() const{
  return _rotMatrix;
}

const Matrix3d & Grain3d::getTempRotMatrix() const{
  return _tempRotMatrix;
}

const Matrix3d & Grain3d::getInertiaTensor() const{
  return _inertiaTensor;
}

const Matrix3d Grain3d::getGlobalI() const{
  return _rotMatrix * _inertiaTensor * _rotMatrix.transpose();
}

const double & Grain3d::getMass() const{
  return _mass;
}

const int & Grain3d::getId() const{
  return _id;
}

const Vector3d & Grain3d::getOmegaGlobal() const{
  return _omegaGlobal;
}

const vector<Vector3d> & Grain3d::getPointList() const {
  return _pointList;
}

const Vector3d & Grain3d::getMomentInertia() const{
  return _momentInertia;
}

const Levelset3d & Grain3d::getLset() const {
  return _lset;
}

void Grain3d::changeMu(const double & mu){
  _mu = mu;
}

void Grain3d::changeQuat(const Vector4d & quat){
  _quat = quat;
}
void Grain3d::changeQuat(const double quat[]){
  for(size_t i = 0; i < 4; ++i) _quat(i) = quat[i];
}

void Grain3d::addVelocityCompression(const Vector3d & vel, const Vector3d & normal, const int & size){
  if(size > 1){
    _velocity += vel - contactDamping*vel.dot(normal)*normal;
  }else{
    _velocity += vel;
  }
}

void Grain3d::addVelocitySeparation(const Vector3d & vel){
  _velocity += vel;
}

void Grain3d::addVelocity(const Vector3d & vel){
  _velocity += vel;
}

void Grain3d::addOmegaGlobal(const Vector3d & omg){
  _omegaGlobal += omg;
  _omega = _rotMatrix.transpose()*_omegaGlobal;
}

void Grain3d::addVelocityAcceleration(const Vector3d & acc, double dt){
  _velocity += acc * dt;
}


void Grain3d::updateQuat(const Matrix3d & R) {
  double tr = R.trace();
  _quat(3) = sqrt(tr+1.)/2.;
  _quat(0) = (R(0,2)-R(2,0))/(4*_quat(3));
  _quat(1) = (R(1,2)-R(2,1))/(4*_quat(3));
  _quat(2) = (R(0,1)-R(1,0))/(4*_quat(3));
}

void Grain3d::updateRotationMatrix(const Vector4d & quat) {
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

void Grain3d::updateTempRotationMatrix(const Vector4d & quat) {
  _tempRotMatrix(0,0) = -quat(0)*quat(0) + quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
  _tempRotMatrix(0,1) = -2*(quat(0)*quat(1) - quat(2)*quat(3));
  _tempRotMatrix(0,2) =  2*(quat(1)*quat(2) + quat(0)*quat(3));
  _tempRotMatrix(1,0) = -2*(quat(0)*quat(1) + quat(2)*quat(3));
  _tempRotMatrix(1,1) =  quat(0)*quat(0) - quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
  _tempRotMatrix(1,2) = -2*(quat(0)*quat(2) - quat(1)*quat(3));
  _tempRotMatrix(2,0) =  2*(quat(1)*quat(2) - quat(0)*quat(3));
  _tempRotMatrix(2,1) = -2*(quat(0)*quat(2) + quat(1)*quat(3));
  _tempRotMatrix(2,2) = -quat(0)*quat(0) - quat(1)*quat(1) + quat(2)*quat(2) + quat(3)*quat(3);
}

void Grain3d::setRealGrain(bool realGrain){
  _realGrain = realGrain;
}

const bool & Grain3d::isRealGrain() const{
  return _realGrain;
}

static std::unique_ptr<Grain3d[]> grainsWorld; //each process have a copy of all grains

struct ContactInfo{
public:
  int      _master;
  int      _slave;
  Vector3d _masterR;
  Vector3d _slaveR;
  Vector3d _normal;
  Vector3d _velocity;
  double   _normalV;
  double   _energy;

  ContactInfo(){
    _master   = -1;
    _slave    = -1;
    _masterR  << 0.,0.,0.;
    _slaveR   << 0.,0.,0.;
    _normal   << 0.,0.,0.;
    _velocity << 0.,0.,0.;
    _normalV  = 0.;
    _energy   = 0.;
  }

  ContactInfo(int master, int slave,
              Vector3d masterR, Vector3d slaveR, Vector3d normal, Vector3d velocity,
              double normalV, double energy):
              _master(master), _slave(slave),
              _masterR(masterR), _slaveR(slaveR), _normal(normal), _velocity(velocity),
              _normalV(normalV), _energy(energy){}

  ContactInfo(const ContactInfo & other){
    _master   = other._master;
    _slave    = other._slave;
    _masterR  = other._masterR;
    _slaveR   = other._slaveR;
    _normal   = other._normal;
    _velocity = other._velocity;
    _normalV  = other._normalV;
    _energy   = other._energy;
  }

  ContactInfo & operator=(const ContactInfo & other){
    _master   = other._master;
    _slave    = other._slave;
    _masterR  = other._masterR;
    _slaveR   = other._slaveR;
    _normal   = other._normal;
    _velocity = other._velocity;
    _normalV  = other._normalV;
    _energy   = other._energy;
    return *this;
  }

  void computeRelVelocity(){
    _velocity = grainsWorld[_master].getVelocity() - grainsWorld[_slave].getVelocity() + grainsWorld[_master].getOmegaGlobal().cross(_masterR) - grainsWorld[_slave].getOmegaGlobal().cross(_slaveR);
    _normalV = _velocity.dot(_normal);
  }

  void computeRelVelocityThreshold(){
    double masterMass = grainsWorld[_master].getMass();
    double slaveMass = grainsWorld[_slave].getMass();
    _velocity = grainsWorld[_master].getVelocity() - grainsWorld[_slave].getVelocity() + grainsWorld[_master].getOmegaGlobal().cross(_masterR) - grainsWorld[_slave].getOmegaGlobal().cross(_slaveR);
    _normalV = _velocity.dot(_normal);
    double t = masterMass < slaveMass / 10. ? 0.1 : threshold;
    if( abs(_normalV) < t ) { _normalV = 0.; }
  }

  void computeRelVelocityWithBalls(const Vector3d & ballVel){
    _velocity = grainsWorld[_master].getVelocity() - ballVel + grainsWorld[_master].getOmegaGlobal().cross(_masterR);
    _normalV = _velocity.dot(_normal);
  }

  void computeRelVelocityWithBallsThreshold(const Vector3d & ballVel){
    _velocity = grainsWorld[_master].getVelocity() - ballVel + grainsWorld[_master].getOmegaGlobal().cross(_masterR);
    _normalV = _velocity.dot(_normal);
    if( abs(_normalV) < threshold ) { _normalV = 0.; }
  }


  void addEnergy(const double & e){
    _energy += e;
  }

  void releaseEnergy(const double & e){
    _energy -= e;
    if(abs(_energy) < threshold) { _energy = 0.; }
  }

  void printContact(){
    cout<<" normal: "<<_normal(0)<<" "<<_normal(1)<<" "<<_normal(2)<<" normal velocity: "<<_normalV<<" master: "<<_master<<" slave: "<<_slave<<endl;
  }

};

#endif /*GRAIN3D_HPP_*/
