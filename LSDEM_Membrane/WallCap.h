/*
 * WallCap.h
 *
 *  Modified on: July 21, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef WALLCAP_H_
#define WALLCAP_H_
#include "definitions.h"
#include "Grain3d.h"

extern int openmpThreads;
extern int dynamicSchedule;

class WallCap{
public:
  WallCap(){
    _position << 0.,0.,0.;
    _normal << 0.,0.,0.;
    _normalGlobal << 0.,0.,0.;
    _center << 0.,0.,0.;
    _ramPos << 0.,0.,0.;
    _d = 0;
    _kn = 0;
    _mu = 0;
    _id = -1;
    _velocity << 0.,0.,0.;
    _omega << 0.,0.,0.,0.;
    _force << 0.,0.,0.;
    _moment << 0.,0.,0.;
    _damping = 0.;
    _height = 0.;
    _radius = 0;
    _density = 0;
    _momentInertia << 0.,0.,0.;
    _dOmega << 0.,0.,0.;
    _rotMatrix << 1.,0.,0.,0.,1.,0.,0.,0.,1.;
    _quat << 0.,0.,0.,1.;
    _omegaGlobal << 0.,0.,0.;
  }

  WallCap(const Vector3d & normal, const Vector3d & position, const double & kn, const double & mu, const int & id, const double & damping, const double & height, const double & radius, const double & density):
      _normal(normal), _position(position), _kn(kn), _mu(mu), _id(id), _damping(damping), _height(height), _radius(radius), _density(density){
    _d = -normal.dot(position); //ax + by + cz + d = 0,
    _velocity << 0.,0.,0.;
    _omega << 0.,0.,0.;
    _force << 0.,0.,0.;
    _moment << 0.,0.,0.;
    _ks = 0.8*_kn;
    _mass = M_PI*_radius*_radius*_height*_density;
    _momentInertia << _mass*(3*_radius*_radius+_height*_height)/12. + _mass*_height*_height/4., _mass*(3*_radius*_radius+_height*_height)/12 + _mass*_height*_height/4., _mass*_radius*_radius/2.;
    _dOmega << 0.,0.,0.;
    _rotMatrix << 1.,0.,0.,0.,1.,0.,0.,0.,1.;
    _quat << 0.,0.,0.,1.;
    _omegaGlobal = _rotMatrix * _omega;
    _normalGlobal = _rotMatrix * _normal;
    _normalGlobal = _normalGlobal/_normalGlobal.norm();
    _center = _position - _height/2.*_normalGlobal;
    _ramPos = _position - _height*_normalGlobal;
  }

  bool bCircleCheck(const std::unique_ptr<Grain3d> & grain) const {
    if(_normalGlobal.dot(grain->getPosition())+_d > grain->getRadius()) { return false; }
    return true;
  }

  void findWallForceMoment(std::unique_ptr<Grain3d> & grain, const double & dt, Vector6d & stressVoigt) {
    Vector3d force(0., 0., 0.);
    Vector3d moment(0., 0., 0.);
    Vector3d capMoment(0., 0., 0.);
    double   penetration;
    Vector3d df;
    Vector3d v;
    Vector3d Fs;
    Vector3d ds;
    double   Fsmag;
    Vector3d ptcm; // point to center of mass (grain)
    Vector3d ptcc; // point to ram
    Vector3d k;
    double sint, cost;
    vector<Vector3d> pointList = grain->getPointList();
    #pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads) private(penetration, df, v, Fs, ds, Fsmag, ptcm, ptcc, k, sint, cost) shared(pointList) reduction(+: force, moment, capMoment)
    for (int i = 0; i < pointList.size(); ++i) {
      penetration = _normalGlobal.dot(pointList[i]) + _d;
      if (penetration < 0 && (pointList[i] - _position).norm() < _radius) {
        ptcm = pointList[i] - grain->getPosition();
        ptcc = pointList[i] -_ramPos;
        df = -penetration*_normalGlobal*_kn;
        v = grain->getVelocity() + grain->getOmegaGlobal().cross(ptcm) - _velocity - _omegaGlobal.cross(ptcc);
        ds = (v-v.dot(_normalGlobal)*_normalGlobal)*dt;
        if(grain->getNormal(i).norm() > 0.){
					k = (-_normalGlobal).cross(grain->getNormal(i));
					sint = k.norm(); cost = sqrt(1-sint*sint);
          if (sint > 0.){ k = k/(sint+std::numeric_limits<double>::epsilon()); }
					if(isnan(cost)){cost = 0.;}
					Vector3d lastShear = grain->getShear(i);
					grain->changeShear(i, lastShear*cost+k.cross(lastShear)*sint+k*k.dot(lastShear)*(1.0-cost));
				}
        if(grain->getNode(i) != -1 && grain->getNode(i) != _id){
          grain->changeNode(i, -1);
          grain->changeShear(i, Vector3d(0.,0.,0.));
          grain->changeNormal(i, Vector3d(0.,0.,0.));
        }
        grain->changeNode(i, _id);
        grain->addShear(i, -ds*_ks);
        grain->changeNormal(i, -_normalGlobal);
        Fsmag = min(df.norm()*_mu, grain->getShear(i).norm());
        if(Fsmag > 0){
          Fs = Fsmag*grain->getShear(i)/grain->getShear(i).norm();
          grain->changeShear(i, Fs);
        }else{ Fs << 0.,0.,0.; }
        Vector3d branchVec = -ptcm;
        Vector3d grainForce = df + Fs;
        force += grainForce; moment += ptcm.cross(grainForce); capMoment += ptcc.cross(-grainForce);
        stressVoigt(0) += grainForce(0)*branchVec(0);
        stressVoigt(1) += grainForce(1)*branchVec(1);
        stressVoigt(2) += grainForce(2)*branchVec(2);
        stressVoigt(3) += 0.5*(grainForce(1)*branchVec(2) + grainForce(2)*branchVec(1));
        stressVoigt(4) += 0.5*(grainForce(2)*branchVec(0) + grainForce(0)*branchVec(2));
        stressVoigt(5) += 0.5*(grainForce(1)*branchVec(0) + grainForce(0)*branchVec(1));
      }else if(grain->getNode(i) == _id){
        grain->changeNode(i, -1);
        grain->changeShear(i, Vector3d(0.,0.,0.));
        grain->changeNormal(i, Vector3d(0.,0.,0.));
      }
    }
    grain->addForce(force);
    grain->addMoment(moment);
    _force  -= force;
    _moment += capMoment;
  }

  void takeTimeStepStrainControl(const Vector3d & disp, double dt){
    Vector3d principleMoment = _rotMatrix.transpose()*_moment; // global to local
    Vector3d omegaN = _omega;
    for (int i = 0; i < 5; i++) {
  		_dOmega(0) = (principleMoment(0) + omegaN(1)*omegaN(2)*(_momentInertia(1)-_momentInertia(2)) - _damping*_momentInertia(0)*omegaN(0) )*dt/_momentInertia(0);
  		_dOmega(1) = (principleMoment(1) + omegaN(2)*omegaN(0)*(_momentInertia(2)-_momentInertia(0)) - _damping*_momentInertia(1)*omegaN(1) )*dt/_momentInertia(1);
  		_dOmega(2) = (principleMoment(2) + omegaN(0)*omegaN(1)*(_momentInertia(0)-_momentInertia(1)) - _damping*_momentInertia(2)*omegaN(2) )*dt/_momentInertia(2);
  		omegaN = _omega + _dOmega/2.;
  	}
    if(isinf(_omega(0)) || isnan(_omega(0)) || isinf(_omega(1)) || isnan(_omega(1)) || isinf(_omega(2)) || isnan(_omega(2))){
      cout << "Issue in WallCap in strain control goes into infinity..." << endl;
      cout << "omega is: "<<_omega(0)<<" "<<_omega(1)<<" "<<_omega(2)<<endl;
      cout << "velocity is: "<<_velocity(0)<<" "<<_velocity(1)<<" "<<_velocity(2)<<endl;
      _dOmega << 0., 0., 0.;
      cout<<"############## END OUTPUT ################"<<endl;
    }
    _omega += _dOmega;
    updateQuatRotationMatrix(dt);
    _ramPos += disp;
    _ramPos << 0., 0., _ramPos(2); // only allow vertical displacement
    _omegaGlobal = _rotMatrix*_omega;
    _normalGlobal = _rotMatrix*_normal;
    _center = _ramPos + _height/2.*_normalGlobal;
    _position = _ramPos + _height*_normalGlobal;
    _d = -_normalGlobal.dot(_position);
    _force << 0., 0., 0.;
    _moment << 0., 0., 0.;
  }

  void takeTimeStepConsolidate(double pressure, double dt){
    _force += pressure*_radius*_radius*M_PI*_normalGlobal;
    _velocity = 1./(1.+_damping*dt/2.)*((1.-_damping*dt/2.)*_velocity+dt*_force/_mass);
    _ramPos += dt*_velocity;
    _ramPos << 0., 0., _ramPos(2); // only allow vertical displacement
    _center = _ramPos + _height/2.*_normalGlobal;
    _position = _ramPos + _height*_normalGlobal;
    _d = -_normalGlobal.dot(_position);
    _force << 0., 0., 0.;
    _moment << 0., 0., 0.;
  }

  const double getArea() const{
    return M_PI*_radius*_radius;
  }

  const Matrix3d & getRotMatrix() const{
    return _rotMatrix;
  }

  const Vector3d & getPosition() const{
    return _position;
  }

  const Vector3d & getCenter() const{
    return _center;
  }

  const Vector4d & getQuat() const{
    return _quat;
  }

  const Vector3d & getForce() const{
    return _force;
  }

  const Vector3d & getMoment() const{
    return _moment;
  }

  Vector3d & getForceNonConst(){
    return _force;
  }

  Vector3d & getMomentNonConst(){
    return _moment;
  }

  Vector3d getNormalGlobal(){
    return _normalGlobal;
  }

  const Vector3d & getRamPos() const{
    return _ramPos;
  }

  double getRadius() const{
    return _radius;
  }

  void changeForce(const Vector3d & force){
    _force = force;
  }

  void changeMoment(const Vector3d & moment){
    _moment = moment;
  }

  void addForce(const Vector3d & dForce){
    _force += dForce;
  }

  void addMoment(const Vector3d & dMoment){
    _moment += dMoment;
  }

  void setNormalGlobal(Vector3d n){
    _normal = n;
  }

  private:
    Vector3d _normal;        // base
    Vector3d _normalGlobal;
    Vector3d _position;      // 'center' at bottom side
    Vector3d _ramPos;
    Vector3d _center;        // 'center' of whole cylinder cap
    double   _kn;            // wall stiffness
    double   _ks;            // shear stiffness
    double   _mu;
    Vector3d _velocity;
    Vector3d _omega;
    Vector3d _dOmega;
    double   _d;             // such that the plane can be written in form ax + by + cz + d = 0 where (a,b,c) = _normal
    int      _id;
    Vector3d _force;
    Vector3d _moment;
    double   _damping;
    double   _height;
    double   _radius;
    double   _density;
    double   _mass;
    Vector3d _momentInertia;
    Matrix3d _rotMatrix;
    Vector4d _quat;
    Vector3d _omegaGlobal;

    void updateQuatRotationMatrix(double dt){
      _quat << -_quat(0), -_quat(1), -_quat(2), _quat(3) ;
      double Bx = dt/4.*_omega(0);
      double By = dt/4.*_omega(1);
      double Bz = dt/4.*_omega(2);
      double c1 =  _quat(0) + Bz*_quat(1) - Bx*_quat(2) - By*_quat(3);
      double c2 =  -Bz*_quat(0) + _quat(1) - By*_quat(2) + Bx*_quat(3);
      double c3 =  Bx*_quat(0) + By*_quat(1) + _quat(2) + Bz*_quat(3);
      double c4 =  By*_quat(0) - Bx*_quat(1) - Bz*_quat(2) + _quat(3);
      double detB = 1 + 2*Bx*Bx + 2*By*By + 2*Bz*Bz + 2*Bx*Bx*By*By + 2*By*By*Bz*Bz + 2*Bx*Bx*Bz*Bz + Bx*Bx*Bx*Bx + By*By*By*By + Bz*Bz*Bz*Bz;
      double bfac = (1 + Bx*Bx + By*By + Bz*Bz)/detB;
      _quat(0) = ( c1 + c2*Bz - c3*Bx - c4*By)*bfac;
      _quat(1) = (-c1*Bz + c2 - c3*By + c4*Bx)*bfac;
      _quat(2) = ( c1*Bx + c2*By + c3 + c4*Bz)*bfac;
      _quat(3) = ( c1*By - c2*Bx - c3*Bz + c4)*bfac;
      _quat = _quat/_quat.norm(); // normalize _quat and use it to update rotation matrix
      _quat << -_quat(0), -_quat(1), -_quat(2), _quat(3); // convert _quaternion back (going from local to global, which is what it usually is)

      _rotMatrix(0,0) = -_quat(0)*_quat(0) + _quat(1)*_quat(1) - _quat(2)*_quat(2) + _quat(3)*_quat(3);
      _rotMatrix(0,1) = -2*(_quat(0)*_quat(1) - _quat(2)*_quat(3));
      _rotMatrix(0,2) =  2*(_quat(1)*_quat(2) + _quat(0)*_quat(3));
      _rotMatrix(1,0) = -2*(_quat(0)*_quat(1) + _quat(2)*_quat(3));
      _rotMatrix(1,1) =  _quat(0)*_quat(0) - _quat(1)*_quat(1) - _quat(2)*_quat(2) + _quat(3)*_quat(3);
      _rotMatrix(1,2) = -2*(_quat(0)*_quat(2) - _quat(1)*_quat(3));
      _rotMatrix(2,0) =  2*(_quat(1)*_quat(2) - _quat(0)*_quat(3));
      _rotMatrix(2,1) = -2*(_quat(0)*_quat(2) + _quat(1)*_quat(3));
      _rotMatrix(2,2) = -_quat(0)*_quat(0) - _quat(1)*_quat(1) + _quat(2)*_quat(2) + _quat(3)*_quat(3);
    }
};
#endif /*WALLCAP_H_*/
