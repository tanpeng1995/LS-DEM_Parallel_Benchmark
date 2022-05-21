/*
 * UTILITIES_H_
 * it is a pity that I made almost all variables global, this is an ugly programming style.
 *
 *  Created on: May 6, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef UTILITIES_H_
#define UTILITIES_H_
#include "Grain3d.h"
#include "WallBalls.h"
#include "WallPlane.h"
#include "WallCap.h"
#include "WallCylinder.h"

//global variables
bool restart;
bool massScaling;
string testName;
string outdirName;
int tmax;
int tConsolidate;
double timeFac;
int outputFrequency;
double scalingFactor;
int num_grains;
double grainDensity;
double kn;
double ks;
double mu;
double gDamping;
double worldWidth;
double worldLength;
double worldHeight;
Vector3d capPosition;
Vector3d capNormalGlobal;
bool dynamic_binning;
int binning_frequency;
double initHeight;
double initRadius;
double ballRadius;
double pressure;
double wallDensity;
double wallToGrainStiffness;
double kmn;
double kms;
double wallFriction;
double dHeight;
double capHeight;
double capRadius;
double capDensity;
double nowCapPressure = 0.;
double prevCapPressure = 0.;
double planePressure = 0.;
double minMass = std::numeric_limits<double>::max();
double maxMass = std::numeric_limits<double>::min();
double maxR = std::numeric_limits<double>::min();
double minX = INFINITY, maxX = -INFINITY, minY = INFINITY, maxY = -INFINITY, minZ = INFINITY, maxZ = -INFINITY;

std::unique_ptr<Grain3d[]> grainsWorld; //each process have a copy of all grains
std::shared_ptr<WallBalls> wallBalls; //specified to extra_proc if exists
std::shared_ptr<WallBalls> externalWall;
std::shared_ptr<WallPlane> wallPlane;
std::shared_ptr<WallCap> wallCap;
std::shared_ptr<WallCylinder> wallCylinder;
std::vector<set<int>> bins;

Vector6d stressVoigt;
#endif
