#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_
#include "definitions.hpp"
//global variables
static bool restart;
static string testName;
static int tmax;
static int tConsolidate;
static double dt;
static double gDamping;
static double threshold;
static int outputFrequency;
static double scalingFactor;
static int num_grains;
static double grainDensity;
static double mu;
static double wallFriction;
static double wallStiffness;
static double dHeight;
static double ballRadius;
static double pressure;
static double ballDensity;
static double kn;
static double kmn;

enum DIRECTIONS {LEFT=0, RIGHT, FRONT, BACK, UP, DOWN};
#endif
