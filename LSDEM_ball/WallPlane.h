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

class WallPlane{
public:
  WallPlane(){ _d = 0; _kn = 0; _mu = 0; _id = -1;}
  WallPlane(const Vector3d & normal, const Vector3d & position, const double & kn, const double & mu, const double & radius, const int & id):
      _normal(normal), _position(position), _kn(kn), _mu(mu), _radius(radius), _id(id){
    _normal = _normal/_normal.norm();
    _d = -(normal.dot(position)); //ax + by + cz + d = 0
    _velocity << 0,0,0; _ks = 0.8*_kn;
    _force << 0.,0.,0.;
  }

  bool bCircleCheck(Grain3d & grain) const {
    if(_normal.dot(grain.getPosition())+_d > grain.getRadius()) return false;
    return true;
  }

  void findWallForceMoment(Grain3d & grain, const double & dt, Vector6d & stressVoigt) {
    Vector3d grainForce(0., 0., 0.);
    Vector3d newNodeNormal(0.,0.,0.);
    Vector3d newNodeShear(0.,0.,0.);
    Vector3d branchVec = grain.getRadius()*_normal;
    int      newNodeContact;
    double   penetration;
    Vector3d v;
    Vector3d Fn(0.,0.,0.);
    Vector3d Fs(0.,0.,0.);
    Vector3d ds;
    double   Fsmag;
    Vector3d k;
    double sint, cost;
    Vector3d contactPos = grain.getPosition() - _normal*grain.getRadius();
    bool found;
    grain.getIfContact(_id, found);
    penetration = _normal.dot(contactPos) + _d;
    if (penetration < 0. && (contactPos - _position).norm() < _radius) {
      /* normal force */
      Fn = -penetration*_normal*_kn;
      /* tangential force */
      v = grain.getVelocity() - _velocity;
      ds = (v-v.dot(_normal)*_normal)*dt;
      if(found){
        Vector3d nodeNormal, nodeShear;
        std::tie(nodeNormal, nodeShear) = grain.getFrictionTuple(_id);
        k = _normal.cross(nodeNormal);
        sint = k.norm(); cost = sqrt(1-sint*sint); if(isnan(cost)){cost = 0.;}
        k = k/(sint+std::numeric_limits<double>::epsilon());
        newNodeShear = nodeShear*cost + k.cross(nodeShear)*sint + k*k.dot(nodeShear)*(1.0-cost);
      }
      newNodeNormal = _normal;
      newNodeShear -= ds*_ks;
      newNodeContact = _id;
      Fsmag = min(Fn.norm()*_mu, newNodeShear.norm());
      if(Fsmag > 0){ newNodeShear = Fsmag*newNodeShear/newNodeShear.norm(); Fs = newNodeShear; }
      grain.changeFrictionTuple(newNodeContact, std::make_tuple(newNodeNormal, newNodeShear));
      grainForce = Fn + Fs;
      stressVoigt(0) += grainForce(0)*branchVec(0);
      stressVoigt(1) += grainForce(1)*branchVec(1);
      stressVoigt(2) += grainForce(2)*branchVec(2);
      stressVoigt(3) += 0.5*(grainForce(1)*branchVec(2) + grainForce(2)*branchVec(1));
      stressVoigt(4) += 0.5*(grainForce(2)*branchVec(0) + grainForce(0)*branchVec(2));
      stressVoigt(5) += 0.5*(grainForce(1)*branchVec(0) + grainForce(0)*branchVec(1));
    }else if(found){
      grain.eraseFrictionTuple(_id);
    }
    grain.addForce(grainForce);
    _force  -= grainForce;
  }

  Vector3d & getForceNonConst(){
    return _force;
  }

  const Vector3d & getForce(){
    return _force;
  }

  Vector3d & getNormalGlobal(){
    return _normal;
  }

  void takeTimeStepStrainControl(const Vector3d & amount, double dt) {
    _position += amount;
    _d = -(_normal.dot(_position));
    _force << 0.,0.,0.;
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
    double   _radius;
    int      _id;
    Vector3d _force;
};
# endif /*WALLPLANE_H_*/
