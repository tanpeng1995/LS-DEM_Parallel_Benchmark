/*
 * WallPlane.h
 *
 *  Created on: August 25, 2014
 *      Author: Reid Kawamoto (Caltech)
 *  Modified on: May 6, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef WALLPLANE_H_
#define WALLPLANE_H_
#include "definitions.h"
#include "Grain3d.h"

extern int openmpThreads;
extern int dynamicSchedule;

class WallPlane{
public:
  WallPlane(){ _d = 0; _kn = 0; _mu = 0; _id = -1;}
  WallPlane(const Vector3d & normal, const Vector3d & position, const double & kn, const double & ks, const double & mu, const int & id):
      _normal(normal), _position(position), _kn(kn), _ks(ks), _mu(mu), _id(id){
    _normal = _normal/_normal.norm();
    _d = -(normal.dot(position)); //ax + by + cz + d = 0
    _velocity << 0,0,0;
  }

  bool bCircleCheck(const std::unique_ptr<Grain3d> & grain) const {
    if(_normal.dot(grain->getPosition())+_d > grain->getRadius() || _normal.dot(grain->getPosition())+_d < 0.) return false;
    return true;
  }

  bool bCircleCheckNonPenetration(const std::unique_ptr<Grain3d> & grain) const {
    if(_normal.dot(grain->getPosition())+_d > grain->getRadius()) return false;
    return true;
  }

  void findWallForceMoment(std::unique_ptr<Grain3d> & grain, const double & dt) {
    Vector3d force(0.,0.,0);
    Vector3d moment(0.,0.,0);
    double   penetration;
    Vector3d df;
    Vector3d v;
    Vector3d Fs;
    Vector3d ds;
    double   Fsmag;
    Vector3d ptcm;
    vector<Vector3d> pointList = grain->getPointList();
    #pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads) private(penetration, df, v, Fs, ds, Fsmag, ptcm) shared(pointList) reduction(+:force, moment)
    for (int i = 0; i < pointList.size(); ++i) {
      penetration = _normal.dot(pointList[i]) + _d;
      if (penetration < 0) {
        ptcm = pointList[i] - grain->getPosition();
  			ptcm = ptcm + penetration*ptcm/ptcm.norm();
        df = -penetration*_normal*_kn;
        v = grain->getVelocity() + grain->getOmegaGlobal().cross(ptcm) - _velocity;
        ds = (v-v.dot(_normal)*_normal)*dt;
        if(grain->getNode(i) != -1 && grain->getNode(i) != _id){
          grain->changeNode(i, -1);
          grain->changeShear(i, Vector3d(0.,0.,0.));
          grain->changeNormal(i, Vector3d(0.,0.,0.));
        }
        grain->changeNode(i, _id);
        grain->addShear(i, -ds*_ks);
        grain->changeNormal(i, -_normal);
        Fsmag = min(df.norm()*_mu, grain->getShear(i).norm());
        if(Fsmag > 0){
          Fs = Fsmag*grain->getShear(i)/grain->getShear(i).norm();
          grain->changeShear(i, Fs);
        }else{ Fs << 0.,0.,0.; }
        Vector3d branchVec = -ptcm;
        Vector3d grainForce = df + Fs;
        force += grainForce; moment += ptcm.cross(grainForce);
      }
      else if(grain->getNode(i) == _id){
        grain->changeNode(i, -1);
        grain->changeShear(i, Vector3d(0.,0.,0.));
        grain->changeNormal(i, Vector3d(0.,0.,0.));
      }
    }
    grain->addForce(force);
    grain->addMoment(moment);
  }

  Vector3d & getNormalGlobal(){
    return _normal;
  }

  void takeTimeStepStrainControl(const Vector3d & amount, double dt) {
    _position += amount;
    _d = -(_normal.dot(_position));
  }

  void takeTimeStepStrainControl(const double & amount) {
    _position(2) += amount;
    _d = -(_normal.dot(_position));
  }

  void changeHeight(const double & newHeight){
    _position << 0.,0.,newHeight;
    _d = -(_normal.dot(_position));
  }

  void increaseHeight(const double & dh){
    _position(2) += dh;
    _d = -(_normal.dot(_position));
  }

  const double & getHeight() const{ return _position(2); }
  private:
    Vector3d _normal;        // outward normal of wall (normalized to magnitude 1 in constructor)
    Vector3d _position;      // 'center' of wall
    double   _kn;            // wall stiffness
    double   _ks;            // shear stiffness
    double   _mu;
    Vector3d _velocity;      // velocity in the normal direction
    double   _d;             // such that the plane can be written in form ax + by + cz + d = 0 where (a,b,c) = _normal
    int      _id;
};
# endif /*WALLPLANE_H_*/
