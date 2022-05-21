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


	bool bCircleCheck(const Grain3d & grain) const {
		Vector3d pos; double gRadius, gRadial;
		pos = grain.getPosition();
		gRadius = grain.getRadius();
		gRadial = sqrt(pos(0)*pos(0) + pos(1)*pos(1));
		if (gRadial+gRadius > _radius) { return true; }
		return false;
	}

	void findWallForceMoment(Grain3d & grain, double dt, Vector6d & stressVoigt) const {
		Vector3d grainPos = grain.getPosition();
		Vector3d normal(grainPos(0), grainPos(1), 0.);
		double   gRadial = sqrt(grainPos(0)*grainPos(0) + grainPos(1)*grainPos(1));
		double   penetration = _radius - gRadial - grain.getRadius();
		Vector3d Fn = penetration*normal*_kn;
		grain.addForce(Fn);
	}

	const double & getHeight() const {
		return _height;
	}
	const double & getRadius() const {
		return _radius;
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
