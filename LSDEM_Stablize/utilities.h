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
#include "WallPlane.h"
#include "WallCylinder.h"
#include "WallBowl.h"

const int LIMIT = 1; //least number of grains per processor
const int SIZE  = 200; //number of grains per bins
const double CUTOFF = 75; // notice here can introduce a bug... SIZE is associated with CUTOFF..

//VERY IMPORTANT: when a struct contains Vector4d for mpi communication, better to use array instead.
//VERY IMPORTANT: when a struct contains Vector4d for mpi communication, better to use array instead.
struct data_t{
  int      _id;
  double   _quat[4];
  Vector3d _position;
  Vector3d _velocity;
  Vector3d _omega;
  data_t(){
    _id = 0; _position << 0.,0.,0.; _velocity << 0.,0.,0.; _omega << 0.,0.,0.;
    for(int i = 0; i < 4; ++i) _quat[i] = 0.;
  }
  data_t(const int& id, const Vector4d& quat, const Vector3d & position, const Vector3d & velocity, const Vector3d & omega):
    _id(id), _position(position), _velocity(velocity), _omega(omega){ for(int i = 0; i < 4; ++i) _quat[i] = quat(i);}
  data_t(const data_t & other){
    _id = other._id; _position = other._position; _velocity = other._velocity; _omega = other._omega;
    for(int i = 0; i < 4; ++i) _quat[i] = other._quat[i];
  }
  ~data_t(){
    _id = -1; _position << 0.,0.,0.; _velocity << 0.,0.,0.; _omega << 0.,0.,0.;
    for(int i = 0; i < 4; ++i) _quat[i] = 0.;
  }
  data_t& operator=(const data_t & other){
    _id = other._id; _position = other._position; _velocity = other._velocity; _omega = other._omega;
    for(int i = 0; i < 4; ++i) _quat[i] = other._quat[i];
    return *this;
  }
};


struct complete_data_t{
  int      _id;
  int      _npoints;
  double   _quat[4];
  Vector3d _position;
  Vector3d _velocity;
  Vector3d _omega;
  complete_data_t(){
    _id       = 0; _npoints  = 0; _position << 0.,0.,0.; _velocity << 0.,0.,0.; _omega    << 0.,0.,0.;
    for(int i = 0; i < 4; ++i) _quat[i] = 0.;
  }

  complete_data_t(const int & id, const int & npoints, const Vector4d& quat,
    const Vector3d & position, const Vector3d & velocity, const Vector3d & omega):
    _id(id), _npoints(npoints), _position(position), _velocity(velocity), _omega(omega){
      for(int i = 0; i < 4; ++i) _quat[i] = quat(i);
    }

  ~complete_data_t(){
    _id = -1; _npoints = -1; _position << 0.,0.,0.; _velocity << 0.,0.,0.; _omega << 0.,0.,0.;
    for(int i = 0; i < 4; ++i) _quat[i] = 0.;
  }

  complete_data_t(const complete_data_t & other){
    _id = other._id; _npoints = other._npoints; _position = other._position; _velocity = other._velocity; _omega = other._omega;
    for(int i = 0; i < 4; ++i) _quat[i] = other._quat[i];
  }

  complete_data_t & operator=(const complete_data_t & other){
    _id = other._id; _npoints = other._npoints; _position = other._position; _velocity = other._velocity; _omega = other._omega;
    for(int i = 0; i < 4; ++i) _quat[i] = other._quat[i];
    return *this;
  }
};

struct output_data{
  int _id;
  double _quat[4];
  Vector3d _position;
  output_data(){
    _id = -1; _position << 0.,0.,0.;
    for(int i = 0; i < 4; ++i) _quat[i] = 0.;
  }

  output_data(const int & id, const Vector4d & quat, const Vector3d & position): _id(id), _position(position){
    for(int i = 0; i < 4; ++i) _quat[i] = quat(i);
  }

  ~output_data(){
    _id = -1; _position << 0.,0.,0.;
    for(int i = 0; i < 4; ++i) _quat[i] = 0.;
  }

  output_data(const output_data & other){
    _id = other._id; _position = other._position;
    for(int i = 0; i < 4; ++i) _quat[i] = other._quat[i];
  }

  output_data& operator=(const output_data & other){
    _id = other._id; _position = other._position;
    for(int i = 0; i < 4; ++i) _quat[i] = other._quat[i];
    return * this;
  }
};

//global variables
string testName;
string InitState;
string outdirName;
bool massScaling;
int tmax;
double timeFac;
int outputFrequency;
bool dynamicBinning;
int  binningFrequency;
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
double initHeight;
double initRadius;
double wallFriction;
int openmpThreads;
int dynamicSchedule;

double totalMass = 0.;

double minMass = std::numeric_limits<double>::max();
double maxMass = std::numeric_limits<double>::min();
double maxR = std::numeric_limits<double>::min();
double minX = INFINITY, maxX = -INFINITY, minY = INFINITY, maxY = -INFINITY, minZ = INFINITY, maxZ = -INFINITY;

std::vector<std::unique_ptr<Grain3d>> grainsWorld;
std::unique_ptr<int[]> grainIdList;
std::unique_ptr<int[]> binIdList;
std::unique_ptr<int[]> belongToThisRank; //0: not in this rank; 1: real particles; 2: ghost particles;
std::unique_ptr<int[]> firstGrainIdList;
std::unique_ptr<WallPlane> bottomPlane;
std::unique_ptr<WallPlane> topPlane;
std::unique_ptr<WallCylinder> wallCylinder;
std::unique_ptr<WallBowl> wallBowl;
std::vector<set<int>> bins;

MPI_Datatype DATA;// MPI_Datatype for data_t
MPI_Datatype COMPLETE_DATA; // for grains leaving host processors
MPI_Datatype OUTPUT_DATA; //for output visualization
MPI_Datatype BIN; // contain fixed number of DATA

int USED; //needed
int dimX, dimY, dimZ; // number of processor in each dimension
int totalBinsInDimX; // total number of bins in X dimension, does not include ghost bins
int totalBinsInDimY; // total number of bins in Y dimension, does not include ghost bins
int totalBinsInDimZ; // total number of bins in Z dimension, does not include ghost bins
int binsInDimX; //number of bins in X dimension for a processor, include 2 ghost bins layer
int binsInDimY; //number of bins in Y dimension for a processor, include 2 ghost bins layer
int binsInDimZ; //number of bins in Z dimension for a processor, include 2 ghost bins layer
int totalBinsPerProc; //total number of bins in a processor

//used to build cartesian grid, this can avoid edge case
enum DIRECTIONS {LEFT, RIGHT, FRONT, BACK, DOWN, UP};
MPI_Comm USING_COMM, CART_COMM;
int neighbor[6];
MPI_Group world_group, new_group, boundary_group;
//used to move grain leaving host processors
//int nHoriz, nVert; //number of bonded spherical walls in horizontal and vertical direction
//vector<Vector3d> ballPositions; //all positions of walls

vector<vector<complete_data_t>> packer;
vector<vector<Vector3d>> packerNodeShears;
vector<vector<int>> packerNodeContact;
vector<vector<Vector3d>> packerNodeNormals;

int * countSendGrain;
int * dispSendGrain;
int * countRecvGrain;
int * dispRecvGrain;

int * countSendVector3d;
int * dispSendVector3d;
int * countRecvVector3d;
int * dispRecvVector3d;

int * countSendInt;
int * dispSendInt;
int * countRecvInt;
int * dispRecvInt;

int totalGrainRecv;
int totalNodeRecv;

vector<complete_data_t> allToAllSendGrains;
vector<complete_data_t> allTOAllRecvGrains;
vector<Vector3d> allToAllSendNodeShears;
vector<Vector3d> allToAllRecvNodeShears;
vector<int> allToAllSendNodeContact;
vector<int> allToAllRecvNodeContact;
vector<Vector3d> allToAllSendNodeNormals;
vector<Vector3d> allToAllRecvNodeNormals;

//used to exchange ghost bins
data_t (* sendToRight)[SIZE]; //(binsInDimY-2)*(binsInDimZ-2)
data_t (* sendToLeft)[SIZE]; //(binsInDimY-2)*(binsInDimZ-2)
data_t (* sendToBack)[SIZE]; //binsInDimX*(binsInDimZ-2)
data_t (* sendToFront)[SIZE]; //binsInDimX*(binsInDimZ-2)
data_t (* sendToDown)[SIZE]; //binsInDimX*binsInDimY
data_t (* sendToUp)[SIZE]; //binsInDimX*binsInDimY
int * grainNumSendToRight; //(binsInDimY-2)*(binsInDimZ-2)
int * grainNumSendToLeft; //(binsInDimY-2)*(binsInDimZ-2)
int * grainNumSendToBack; //binsInDimX*(binsInDimZ-2)
int * grainNumSendToFront; //binsInDimX*(binsInDimZ-2)
int * grainNumSendToDown; //binsInDimX*binsInDimY
int * grainNumSendToUp; //binsInDimX*binsInDimY
data_t (* recvFromLeft)[SIZE]; //(binsInDimY-2)*(binsInDimZ-2)
data_t (* recvFromRight)[SIZE]; //(binsInDimY-2)*(binsInDimZ-2)
data_t (* recvFromFront)[SIZE]; //binsInDimX*(binsInDimZ-2)
data_t (* recvFromBack)[SIZE]; //binsInDimX*(binsInDimZ-2)
data_t (* recvFromUp)[SIZE]; //binsInDimX*binsInDimY
data_t (* recvFromDown)[SIZE]; //binsInDimX*binsInDimY
int * grainNumRecvFromRight; //(binsInDimY-2)*(binsInDimZ-2)
int * grainNumRecvFromLeft; //(binsInDimY-2)*(binsInDimZ-2)
int * grainNumRecvFromBack; //binsInDimX*(binsInDimZ-2)
int * grainNumRecvFromFront; //binsInDimX*(binsInDimZ-2)
int * grainNumRecvFromDown; //binsInDimX*binsInDimY
int * grainNumRecvFromUp; //binsInDimX*binsInDimY

Vector6d stressVoigt;
#endif
