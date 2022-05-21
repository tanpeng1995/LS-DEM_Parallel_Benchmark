#ifndef WALLPLANE_HPP_
#define WALLPLANE_HPP_
#include "definitions.hpp"
#include "Grain3d.hpp"
#include "utilities.hpp"

#pragma omp declare reduction(+ : Vector3d : omp_out += omp_in ) initializer (omp_priv=Vector3d(0.,0.,0.))

const bool planeDebug = false;

class WallPlane{
public:
  WallPlane();
  WallPlane(const Vector3d & normal, const Vector3d & position, const double & mu, const double & v, const double & kn);
  bool bCircleCheck(const std::unique_ptr<Grain3d> & grain) const;
  bool bTempCircleCheck(const std::unique_ptr<Grain3d> & grain) const;
  void findGrainWallContactImpulse(std::unique_ptr<Grain3d> & grain, double dt, vector<ContactInfo> & contactList);
  void findGrainWallContact(std::unique_ptr<Grain3d> & grain, double dt);
  void movePlaneVelocityControl(double dt);
  void movePlaneStrainControl(const double & amount);
  Vector3d getNormal() const;

  private:
    Vector3d _normal;        // outward normal of wall (normalized to magnitude 1 in constructor)
    Vector3d _position;      // 'center' of wall
    double   _mu;
    Vector3d _velocity;      // velocity in the normal direction
    double   _d;             // such that the plane can be written in form ax + by + cz + d = 0 where (a,b,c) = _normal
    double  _kn;
};

WallPlane::WallPlane(){
  _d = 0;
  _mu = 0;
}

WallPlane::WallPlane(const Vector3d & normal, const Vector3d & position, const double & mu, const double & v, const double & kn):
  _normal(normal), _position(position), _mu(mu), _kn(kn){
  _normal = _normal/_normal.norm();
  _d = -(normal.dot(position)); //ax + by + cz + d = 0
  _velocity = _normal * v;
}

bool WallPlane::bCircleCheck(const std::unique_ptr<Grain3d> & grain) const {
  if(_normal.dot(grain->getPosition())+_d > grain->getRadius()) return false;
  return true;
}

bool WallPlane::bTempCircleCheck(const std::unique_ptr<Grain3d> & grain) const {
  if(_normal.dot(grain->getTempPosition())+_d > grain->getRadius()) return false;
  return true;
}

void WallPlane::findGrainWallContactImpulse(std::unique_ptr<Grain3d> & grain, double dt, vector<ContactInfo> & contactList){
  Vector3d nodeForce(0.,0.,0);
  Vector3d grainMoment(0.,0.,0);
  Vector3d ptcm, v;
  Vector3d netV(0.,0.,0.);
  Vector3d netBV(0.,0.,0.);
  int      nContact = 0;
  double   penetration;
  vector<Vector3d> pointList = grain->getPointList();
  for (int i = 0; i < pointList.size(); ++i) {
    penetration = _normal.dot(pointList[i]) + _d;
    if (penetration < 0) {
      Vector3d ptcm = pointList[i] - grain->getTempPosition() - penetration * _normal;
      v = grain->getVelocity() + grain->getOmegaGlobal().cross(ptcm) - _velocity;
      netBV += ptcm; netV += v; nContact++;
    }
  }

  if(nContact > 0){
    netBV /= nContact; netV /= nContact;
    double relVel = _normal.dot(netV);
    if(relVel < -threshold){
      contactList.emplace_back(grain->getId(), -1, netBV, Vector3d(0.,0.,0.), _normal, netV, relVel, 0.);
    }
  }
}


void WallPlane::findGrainWallContact(std::unique_ptr<Grain3d> & grain, double dt) {
    Vector3d nodeForce(0.,0.,0);
    Vector3d grainMoment(0.,0.,0);
    double   penetration;
    double   minPenetration = std::numeric_limits<double>::max();
    vector<Vector3d> pointList = grain->getPointList();
    for (int i = 0; i < pointList.size(); ++i) {
      penetration = _normal.dot(pointList[i]) + _d;
      minPenetration = minPenetration < penetration ? minPenetration : penetration;
      if (penetration < 0) {
        Vector3d ptcm = pointList[i] - grain->getTempPosition() - penetration * _normal;
        nodeForce = -penetration*_normal*_kn;
        grainMoment += ptcm.cross(nodeForce);
      }
    }

    if(minPenetration < 0.){
      Vector3d dVel = -minPenetration / dt * _normal;
      grain->addVelocity(dVel);
    }

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
      //if(dOmega.norm() > 1. ) dOmega = dOmega / dOmega.norm();
      /* this is important at the first step, initial overlaps */
      omega += dOmega;
      grain->changeOmega(omega);
    }
  }

void WallPlane::movePlaneStrainControl(const double & amount) {
  _position += amount * _normal;
  _d = -(_normal.dot(_position));
}

void WallPlane::movePlaneVelocityControl(double dt)  {
  _position += _velocity * dt;
  _d = -(_normal.dot(_position));
}

Vector3d WallPlane::getNormal() const{
  return _normal;
}

static std::vector<std::unique_ptr<WallPlane>> wallPlanes;

#endif /*WALLPLANE_HPP_*/
