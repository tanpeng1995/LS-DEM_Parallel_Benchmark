#ifndef INIT_SIMULATION_H_
#define INIT_SIMULATION_H_
#include "Grain3d.hpp"
#include "WallPlane.hpp"
#include "WallBalls.hpp"
#include "WallTopography.hpp"
#include "readParameters.hpp"

const int whichWallBall(int slaveId){
  if(slaveId >= num_grains + num_balls) return 1;
  else if(slaveId >= num_grains) return 0;
}

void init_simulation(string indir, FILE * wallfile1, FILE * wallfile2){
  vector<Vector3d> inputPos = readPositionFile(indir + "InitState/positions_"+testName+".dat", num_grains);
  vector<Vector4d> inputRot = readQuaternionFile(indir + "InitState/rotations_"+testName+".dat", num_grains);
  grainsWorld.resize(num_grains);
  wallPlanes.resize(6);

  wallPlanes[LEFT]  = std::make_unique<WallPlane>(Vector3d(1.,0.,0.), Vector3d(0.,0.,0.), wallFriction, 0., wallStiffness);
  wallPlanes[RIGHT] = std::make_unique<WallPlane>(Vector3d(-1.,0.,0.), Vector3d(4000.,0.,0.), wallFriction, 0., wallStiffness);
  wallPlanes[FRONT] = std::make_unique<WallPlane>(Vector3d(0.,-1.,0.), Vector3d(0.,1280.,0.), wallFriction, 0., wallStiffness);
  wallPlanes[BACK]  = std::make_unique<WallPlane>(Vector3d(0.,1.,0.), Vector3d(0.,0.,0.), wallFriction, 0., wallStiffness);
  wallPlanes[UP]    = std::make_unique<WallPlane>(Vector3d(0.,0.,-1.), Vector3d(0.,0.,2400.), wallFriction, 0., wallStiffness);
  wallPlanes[DOWN]  = std::make_unique<WallPlane>(Vector3d(0.,0.,1.), Vector3d(0.,0.,0.), wallFriction, 0., wallStiffness); //tilt 3 degree

  cout<<"finish buiding wallPlane and wallCap..."<<endl;
  wallBalls[0] = std::make_unique<WallBalls>(300., 100., ballRadius, ballDensity, kn, kmn, 0.2247, Vector3d(1200., 100., 200.), num_grains);
  num_balls = wallBalls[0]->getNumWalls();
  wallBalls[1] = std::make_unique<WallBalls>(400., 100., ballRadius, ballDensity, kn, kmn, 0., Vector3d(1300., 50., 200.), num_grains+num_balls);
  wallTopography = 
std::make_unique<WallTopography>("/global/scratch/users/tanpeng/fastmarching/very_big_topography.dat", 0.2247, Vector3d(0.,0.,800.), wallFriction, wallStiffness);
  cout<<"finish buiding wallTopography..."<<endl;

  fprintf(wallfile1, "%d\n", wallBalls[0]->getNumWalls());
  fprintf(wallfile1, "%.4f\n", wallBalls[0]->getBallRadius());
  fflush(wallfile1);

  fprintf(wallfile2, "%d\n", wallBalls[1]->getNumWalls());
  fprintf(wallfile2, "%.4f\n", wallBalls[1]->getBallRadius());
  fflush(wallfile2);

  for(int i = 0; i < num_grains; ++i){
    if(i % 100 == 0) cout<<"finish loading "<<i<<" grains"<<endl;
    stringstream morphfile;
    morphfile << indir << "Morphologies/morph_" << i+1 << ".dat";
    grainsWorld[i] = generateGrainFromFile(inputPos[i], inputRot[i], morphfile.str(), grainDensity, i);
  }
}
#endif /*INIT_SIMULATION_HPP_*/
