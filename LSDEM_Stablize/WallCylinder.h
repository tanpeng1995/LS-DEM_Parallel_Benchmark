/*
 * WallCylinder.h
 *
 *  Created on: Oct 24, 2014
 *      Author: Reid (Caltech)
 *  Modified on: July 10, 2020
 *      Author: Peng TAN (Berkeley)
 */

#ifndef WALLCYLINDER_H_
#define WALLCYLINDER_H_

#include "definitions.h"
#include "Grain3d.h"

extern int openmpThreads;
extern int dynamicSchedule;

class WallCylinder {

public:
	WallCylinder() { _height = 0; _radius = 0; _kn = 0; _ks = 0; _mu = 0; _id = 0;}
	WallCylinder(const double & height, const double & radius, const double & kn, const double & ks, const double & mu, const int & id):
	_height(height), _radius(radius), _kn(kn), _ks(ks), _mu(mu), _id(id) {}


	bool bCircleCheck(const std::unique_ptr<Grain3d> & grain) const {
		Vector3d pos; double gRadius, gRadial;
		pos = grain->getPosition();
		gRadius = grain->getRadius();
		gRadial = sqrt(pos(0)*pos(0) + pos(1)*pos(1));
		if (gRadial+gRadius > _radius) { return true; }
		return false;
	}

	void findWallForceMoment(std::unique_ptr<Grain3d> & grain, double dt) const {
		Vector3d force(0.,0.,0.);
		Vector3d moment(0.,0.,0.);
		double   gRadial;	// radial component of grain's x-y coordinates
		Vector3d pos;		// position of grain points
		double   penetration; // penetration amount into wall
		Vector3d ptcm;		// location of point wrt the grain center of mass
		Vector3d df;	// incremental force on the grain

		Vector3d normal;
		Vector3d v;
		Vector3d ds;
		Vector3d Fs;
		double Fsmag;

		Vector3d k;
		double sint, cost;
		vector<Vector3d> pointList = grain->getPointList();
		#pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads) private(gRadial, pos, penetration, ptcm, df, normal, v, ds, Fs, Fsmag, k, sint, cost) shared(force, moment, grain)
    for (int i = 0; i < pointList.size(); i++) {
      pos = pointList[i];
      gRadial = sqrt(pos(0)*pos(0) + pos(1)*pos(1));
      penetration = _radius - gRadial;
      if (penetration < 0.) {
        ptcm = pos - grain->getPosition();
        normal << pos(0)/gRadial, pos(1)/gRadial, 0;
        df = penetration*normal*_kn;
				force += df; moment += ptcm.cross(df);

				v = grain->getVelocity() + grain->getOmegaGlobal().cross(ptcm);
				ds = (v - v.dot(normal)*normal)*dt;

				if(grain->getNormal(i).norm()>0){
					k = normal.cross(grain->getNormal(i));
					sint = k.norm(); cost = sqrt(1-sint*sint); k = k/(sint+std::numeric_limits<double>::epsilon()); //DBL_MIN is not defined?
					if(isnan(cost)){cost = 0.;}
					Vector3d lastShear = grain->getShear(i);
					grain->changeShear(i, lastShear*cost+k.cross(lastShear)*sint+k*k.dot(lastShear)*(1.0-cost));
				}

				grain->changeNode(i, _id);
        grain->changeShear(i, -ds*_ks);
        grain->changeNormal(i, normal);
        Fsmag = min(df.norm()*_mu, grain->getShear(i).norm());
        if(Fsmag > 0){
          Fs = Fsmag*grain->getShear(i)/grain->getShear(i).norm();
          grain->changeShear(i, Fs);
					force += Fs;
					moment += ptcm.cross(Fs);
        }
      }
    }
		grain->addForce(force);
		grain->addMoment(moment);
	}

	void findWallForceMomentNonFriction(std::unique_ptr<Grain3d> & grain, double dt) const {
		Vector3d force(0.,0.,0.);
		Vector3d moment(0.,0.,0.);
		double   gRadial;	// radial component of grain's x-y coordinates
		Vector3d pos;		// position of grain points
		double   penetration; // penetration amount into wall
		Vector3d ptcm;		// location of point wrt the grain center of mass
		Vector3d df;	// incremental force on the grain

		Vector3d normal;
		Vector3d v;
		Vector3d ds;
		Vector3d Fs;
		double Fsmag;

		Vector3d k;
		double sint, cost;
		vector<Vector3d> pointList = grain->getPointList();
		#pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads) private(gRadial, pos, penetration, ptcm, df, normal, v, ds, Fs, Fsmag, k, sint, cost) shared(force, moment, grain)
    for (int i = 0; i < pointList.size(); i++) {
      pos = pointList[i];
      gRadial = sqrt(pos(0)*pos(0) + pos(1)*pos(1));
      penetration = _radius - gRadial;
      if (penetration < 0.) {
        ptcm = pos - grain->getPosition();
        normal << pos(0)/gRadial, pos(1)/gRadial, 0;
        df = penetration*normal*_kn;
				force += df; moment += ptcm.cross(df);
      }
    }
		grain->addForce(force);
		grain->addMoment(moment);
	}

	const double & getHeight() const {
		return _height;
	}
	const double & getRadius() const {
		return _radius;
	}

	void changeRadius(const double & dr){
		_radius += dr;
	}


private:
	double _height;	// height of wall
	double _radius;	// radius of wall
	double _kn;			// normal stiffness of wall
	double _ks;
	double _mu;			// friction of wall
	int _id;

};



#endif /* WALLCYLINDER_H_ */
