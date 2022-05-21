#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_
#include "definitions.hpp"
#include "utilities.hpp"

const int LIMIT = 1; //least number of grains per processor
const int SIZE  = 200; //number of grains per bins
const double CUTOFF = 50; // notice here can introduce a bug... SIZE is associated with CUTOFF..

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

struct vel_data{
public:
  int _id;
  int _rank;
  Vector3d _velocity;
  Vector3d _omegaGlobal;

  vel_data(): _id(-1), _rank(-1), _velocity(Vector3d(0.,0.,0.)), _omegaGlobal(Vector3d(0.,0.,0.)){}
  vel_data(const int & id, const int & rank, const Vector3d & velocity, const Vector3d & omegaGlobal):
    _id(id), _rank(rank), _velocity(velocity), _omegaGlobal(omegaGlobal){}
  vel_data(const vel_data & other):
    _id(other._id), _rank(other._rank), _velocity(other._velocity), _omegaGlobal(other._omegaGlobal){}
  vel_data & operator=(const vel_data & other){
    _id = other._id;
    _rank = other._rank;
    _velocity = other._velocity;
    _omegaGlobal = other._omegaGlobal;
    return *this;
  }
};

//global variables
static bool restart;
static bool massScaling;
static string testName;
static int tmax;
static double dt;
static double gDamping;
static double contactDamping;
static double threshold;
static double similarityThreshold;
static int outputFrequency;
static double scalingFactor;
static int num_grains;
static int target_grains;
static double grainDensity;
static double mu;
static double worldWidth;
static double worldLength;
static double worldHeight;
static double targetWidth;
static double targetLength;
static double targetHeight;
static double wallStiffness;
static bool dynamic_binning;
static int binning_frequency;
static int openmpThreads;
static int dynamicSchedule;
static double minX = INFINITY, maxX = -INFINITY, minY = INFINITY, maxY = -INFINITY, minZ = INFINITY, maxZ = -INFINITY;
static double minMass, maxMass;

static std::unique_ptr<int[]> grainIdList;
static std::unique_ptr<int[]> binIdList;
static std::unique_ptr<int[]> belongToThisRank; //0: not in this rank; 1: real particles; 2: ghost particles;
static std::unique_ptr<int[]> firstGrainIdList;
static std::vector<set<int>> bins;

static MPI_Datatype DATA;// MPI_Datatype for data_t
static MPI_Datatype OUTPUT_DATA; //for output visualization
static MPI_Datatype BIN; // contain fixed number of DATA
static MPI_Datatype VELOCITY;
static MPI_Datatype CONTACT;

static int USED; //needed
static int dimX, dimY, dimZ; // number of processor in each dimension
static int totalBinsInDimX; // total number of bins in X dimension, does not include ghost bins
static int totalBinsInDimY; // total number of bins in Y dimension, does not include ghost bins
static int totalBinsInDimZ; // total number of bins in Z dimension, does not include ghost bins
static int binsInDimX; //number of bins in X dimension for a processor, include 2 ghost bins layer
static int binsInDimY; //number of bins in Y dimension for a processor, include 2 ghost bins layer
static int binsInDimZ; //number of bins in Z dimension for a processor, include 2 ghost bins layer
static int totalBinsPerProc; //total number of bins in a processor

//used to build cartesian grid, this can avoid edge case
enum DIRECTIONS {LEFT=0, RIGHT, FRONT, BACK, DOWN, UP};
static MPI_Comm CART_COMM;
static int neighbor[6];
static MPI_Group world_group, new_group;
static MPI_Comm USING_COMM;
//used to move grain leaving host processors
//int nHoriz, nVert; //number of bonded spherical walls in horizontal and vertical direction
//vector<Vector3d> ballPositions; //all positions of walls
static std::vector<data_t> grain_to_be_send[6];
static std::vector<data_t> grain_to_be_recv[6];
static int grain_send_count[6];
static int grain_recv_count[6];

//used to exchange ghost bins
static data_t (* sendToRight)[SIZE]; //(binsInDimY-2)*(binsInDimZ-2)
static data_t (* sendToLeft)[SIZE]; //(binsInDimY-2)*(binsInDimZ-2)
static data_t (* sendToBack)[SIZE]; //binsInDimX*(binsInDimZ-2)
static data_t (* sendToFront)[SIZE]; //binsInDimX*(binsInDimZ-2)
static data_t (* sendToDown)[SIZE]; //binsInDimX*binsInDimY
static data_t (* sendToUp)[SIZE]; //binsInDimX*binsInDimY
static int * grainNumSendToRight; //(binsInDimY-2)*(binsInDimZ-2)
static int * grainNumSendToLeft; //(binsInDimY-2)*(binsInDimZ-2)
static int * grainNumSendToBack; //binsInDimX*(binsInDimZ-2)
static int * grainNumSendToFront; //binsInDimX*(binsInDimZ-2)
static int * grainNumSendToDown; //binsInDimX*binsInDimY
static int * grainNumSendToUp; //binsInDimX*binsInDimY
static data_t (* recvFromLeft)[SIZE]; //(binsInDimY-2)*(binsInDimZ-2)
static data_t (* recvFromRight)[SIZE]; //(binsInDimY-2)*(binsInDimZ-2)
static data_t (* recvFromFront)[SIZE]; //binsInDimX*(binsInDimZ-2)
static data_t (* recvFromBack)[SIZE]; //binsInDimX*(binsInDimZ-2)
static data_t (* recvFromUp)[SIZE]; //binsInDimX*binsInDimY
static data_t (* recvFromDown)[SIZE]; //binsInDimX*binsInDimY
static int * grainNumRecvFromRight; //(binsInDimY-2)*(binsInDimZ-2)
static int * grainNumRecvFromLeft; //(binsInDimY-2)*(binsInDimZ-2)
static int * grainNumRecvFromBack; //binsInDimX*(binsInDimZ-2)
static int * grainNumRecvFromFront; //binsInDimX*(binsInDimZ-2)
static int * grainNumRecvFromDown; //binsInDimX*binsInDimY
static int * grainNumRecvFromUp; //binsInDimX*binsInDimY

#endif
