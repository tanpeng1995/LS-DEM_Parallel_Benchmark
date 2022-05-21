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
  WallPlane(const Vector3d & normal, const Vector3d & position, const double & kn, const double & mu, const int & id):
      _normal(normal), _position(position), _kn(kn), _mu(mu), _id(id){
    _normal = _normal/_normal.norm();
    _d = -(normal.dot(position)); //ax + by + cz + d = 0
    _velocity << 0.,0.,0.;
    _ks = 0.8*_kn;
  }

  bool bCircleCheck(Grain3d & grain) const {
    if(_normal.dot(grain.getPosition())+_d > grain.getRadius()) return false;
    return true;
  }

  void findWallForceMoment(Grain3d & grain, const double & dt) {
    Vector3d force(0.,0.,0);
    Vector3d moment(0.,0.,0);
    double   penetration;
    Vector3d df;
    Vector3d v;
    Vector3d Fs;
    Vector3d ds;
    double   Fsmag;
    Vector3d ptcm;
    for (int i = 0; i < grain.getPointList().size(); ++i) {
      penetration = _normal.dot(grain.getPointList()[i]) + _d;
      if (penetration < 0.) {
        ptcm = grain.getPointList()[i] - grain.getPosition();
  			ptcm = ptcm + penetration*ptcm/ptcm.norm();
        df = -penetration*_normal*_kn;
        force += df;
        moment += ptcm.cross(df);
        v = grain.getVelocity() + grain.getOmegaGlobal().cross(ptcm) - _velocity;
        ds = (v-v.dot(_normal)*_normal)*dt;
        grain.getNodeContactNonConst()[i] = _id;
        grain.getNodeShearsNonConst()[i] -= ds*_ks;
        Fsmag = min(df.norm()*_mu, grain.getNodeShears()[i].norm());
        if(Fsmag > 0.){
          Fs = Fsmag*grain.getNodeShears()[i]/grain.getNodeShears()[i].norm();
          grain.getNodeShearsNonConst()[i] = Fs;
          force  += Fs;
          moment += ptcm.cross(Fs);
        }
      }
    }
    grain.addForce(force);
    grain.addMoment(moment);
  }


  void moveHeightStrainControl(const Vector3d & amount) {
    _position += amount;
    _d = -(_normal.dot(_position));
  }

  const double & getHeight() const{
    return _position(2);
  }
  
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
