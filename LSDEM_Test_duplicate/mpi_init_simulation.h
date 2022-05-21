/*
 * MPI_INIT_SIMULATION_H_
 * Created on: June 20, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_INIT_SIMULATION_H_
#define MPI_INIT_SIMULATION_H_
#include "mpi_helper_function.h"
#include "mpi_grains.h"
void init_simulation(int rank, int num_procs, string indir){
  //num_procs-1 because need to leave one processor doint boundary staff
  determine_proc(num_grains, num_procs, USED, dimX, dimY, dimZ);
  sort(dimX, dimY, dimZ); //so that dimZ >= dimX >= dimY
  if(rank == 0){
    cout<<"dimX: "<<dimX<<" dimY: "<<dimY<<" dimZ: "<<dimZ<<endl;
  }
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);

  int ranks[USED];
  for(int i = 0; i < USED; ++i) ranks[i] = i;
  MPI_Group_incl(world_group, USED, ranks, &new_group);
  MPI_Comm_create(MPI_COMM_WORLD, new_group, &USING_COMM);

  //commit MPI_Datatype to exchange data
  const int items1 = 5;
  int lengths1[5] = {1, 4, 3, 3, 3};
  MPI_Datatype types1[5] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offsets1[5];
  offsets1[0] = offsetof(data_t, _id);
  offsets1[1] = offsetof(data_t, _quat);
  offsets1[2] = offsetof(data_t, _position);
  offsets1[3] = offsetof(data_t, _velocity);
  offsets1[4] = offsetof(data_t, _omega);
  MPI_Type_create_struct(items1, lengths1, offsets1, types1, &DATA);
  MPI_Type_commit(&DATA);

  const int items2 = 6;
  int lengths2[6] = {1, 1, 4, 3, 3, 3};
  MPI_Datatype types2[6] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offsets2[6];
  offsets2[0] = offsetof(complete_data_t, _id);
  offsets2[1] = offsetof(complete_data_t, _npoints);
  offsets2[2] = offsetof(complete_data_t, _quat);
  offsets2[3] = offsetof(complete_data_t, _position);
  offsets2[4] = offsetof(complete_data_t, _velocity);
  offsets2[5] = offsetof(complete_data_t, _omega);
  MPI_Type_create_struct(items2, lengths2, offsets2, types2, &COMPLETE_DATA);
  MPI_Type_commit(&COMPLETE_DATA);

  const int items3 = 3;
  int lengths3[3] = {1, 4, 3};
  MPI_Datatype types3[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offsets3[3];
  offsets3[0] = offsetof(output_data, _id);
  offsets3[1] = offsetof(output_data, _quat);
  offsets3[2] = offsetof(output_data, _position);
  MPI_Type_create_struct(items3, lengths3, offsets3, types3, &OUTPUT_DATA);
  MPI_Type_commit(&OUTPUT_DATA);

  MPI_Type_contiguous(SIZE, DATA, &BIN);
  MPI_Type_commit(&BIN);

  if(rank < USED){
    vector<Vector3d> inputPos;
    vector<Vector4d> inputRot;

    inputPos = readPositionFile(indir + "InitState/positions_"+testName+".dat", num_grains);
    inputRot = readQuaternionFile(indir + "InitState/rotations_"+testName+".dat", num_grains);

    grainsWorld = std::make_unique<Grain3d[]>(num_grains);
    grainIdList = std::make_unique<int[]>(num_grains);
    binIdList = std::make_unique<int[]>(num_grains);
    belongToThisRank = std::make_unique<int[]>(num_grains);
    wallPlanes = std::make_unique<WallPlane[]>(6);
    wallPlanes[0] = WallPlane(Vector3d(0.,0.,1.), Vector3d(0.,0.,0.), kn, wallFriction, num_grains); // bottom wall
    wallPlanes[1] = WallPlane(Vector3d(0.,0.,-1.), Vector3d(0.,0.,worldHeight), kn, wallFriction, num_grains+1); // top wall
    wallPlanes[2] = WallPlane(Vector3d(1.,0.,0.), Vector3d(0.,0.,0.), kn, wallFriction, num_grains+2); // right wall
    wallPlanes[3] = WallPlane(Vector3d(-1.,0.,0.), Vector3d(worldLength,0.,0.), kn, wallFriction, num_grains+3); // left wall
    wallPlanes[4] = WallPlane(Vector3d(0.,1.,0.), Vector3d(0.,0.,0.), kn, wallFriction, num_grains+4); // back wall
    wallPlanes[5] = WallPlane(Vector3d(0.,-1.,0.), Vector3d(0.,worldWidth,0.),  kn, wallFriction, num_grains+5); // front wall

    // if rank is normal processor
    int dims[3] = {dimZ, dimY, dimX};
    MPI_Dims_create(USED, 3, dims);
    int periods[3] = {false, false, false};
    int reorder = false;
    MPI_Cart_create(USING_COMM, 3, dims, periods, reorder, &CART_COMM);
    //if this line put outside the loop, will report error, because USING_COMM only includes rank < USED.
    MPI_Cart_shift(CART_COMM, 0, 1, &neighbor[DOWN],  &neighbor[UP]);
    MPI_Cart_shift(CART_COMM, 1, 1, &neighbor[FRONT], &neighbor[BACK]);
    MPI_Cart_shift(CART_COMM, 2, 1, &neighbor[LEFT],  &neighbor[RIGHT]);

    maxX = worldWidth; maxY = worldLength; maxZ = worldHeight;
    minX = 0.; minY = 0.; minZ = 0.;

    worldWidth  = maxX-minX;
    worldLength = maxY-minY;
    worldHeight = maxZ-minZ;

    //Initialize variable
    totalBinsInDimX = ceil(worldWidth/CUTOFF); // total number of bins in X dimension
    totalBinsInDimY = ceil(worldLength/CUTOFF); // total number of bins in Y dimension
    totalBinsInDimZ = ceil(worldHeight/CUTOFF); // total number of bins in Z dimension
    binsInDimX = ceil((double)totalBinsInDimX/dimX) + 2; //number of bins in X dimension for a processor, including 2 ghost
    binsInDimY = ceil((double)totalBinsInDimY/dimY) + 2; //number of bins in Y dimension for a processor, including 2 ghost
    binsInDimZ = ceil((double)totalBinsInDimZ/dimZ) + 2; //number of bins in Z dimension for a processor, including 2 ghost
    totalBinsPerProc = binsInDimX*binsInDimY*binsInDimZ; //number of bins in a processor
    initMPIBuff_static(num_grains, USED);
    initMPIBuff_dynamic(num_grains, SIZE, totalBinsPerProc, binsInDimX, binsInDimY, binsInDimZ);

    //assign particles
    int _idx, _idy, _idz, _bin_id, _proc_id;
    double maxR = 0.;
    int Nx = worldWidth / targetWidth;
    int Ny = worldLength / targetLength;
    int Nz = worldHeight / targetHeight;
    for(int iz = 0; iz < Nz; ++iz){
      for(int iy = 0; iy < Ny; ++iy){
        for(int ix = 0; ix < Nx; ++ix){
          int disp_i = (Nx*Ny*iz + Nx*iy + ix)*target_grains;

          /* BEGAIN DUPLICATING */
          for(int i = 0; i < target_grains; ++i){
            int gid = i + disp_i;
            if(rank == 0 && gid % 100 == 0) cout<<"finish loading "<<gid<<" grains"<<endl;
            find_pos(inputPos[gid], _idx, _idy, _idz);
            _proc_id = to_proc_id(_idx, _idy, _idz);
            _bin_id = to_bin_id(_idx, _idy, _idz);
            stringstream morphfile;
            morphfile << indir << "Morphologies/morph_" << i+1 << ".dat";
            if(_proc_id == rank){
              grainsWorld[gid] = generateGrainFromFile(inputPos[gid], inputRot[gid], morphfile.str(), gid);
              grainsWorld[gid].changeDensity(grainDensity);
              grainsWorld[gid].setRealGrain(true);
              bins[_bin_id].insert(gid);
              binIdList[gid] = _bin_id;

              if(firstGrainIdList[_bin_id] == -1){
                firstGrainIdList[_bin_id] = gid;
              }else{
                int j = firstGrainIdList[_bin_id];
                while(grainIdList[j] != -1) j = grainIdList[j];
                grainIdList[j] = gid;
              }
              belongToThisRank[gid] = 1;
              minMass = min(minMass, grainsWorld[gid].getMass());
              double grainR = grainsWorld[gid].getRadius();
              maxR = maxR > grainR ? maxR : grainR;
            }else if(ghost_bin(_idx, _idy, _idz, rank)){
              grainsWorld[gid] = generateGrainFromFile(inputPos[gid], inputRot[gid], morphfile.str(), gid);
              grainsWorld[gid].changeDensity(grainDensity);
              grainsWorld[gid].setRealGrain(true);
              belongToThisRank[gid] = 2;
            }else{
              grainsWorld[gid] = Grain3d(inputPos[gid], inputRot[gid], gid);
              grainsWorld[gid].setRealGrain(false);
              belongToThisRank[gid] = 0;
            }
          }
          /* END DUPLICATING */

        }
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &minMass, 1, MPI_DOUBLE, MPI_MIN, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &maxR, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    if(rank == 0){ cout<<"largest possible radius is: "<<maxR<<" min mass is: "<<minMass<<endl; }
  }
}
#endif /*MPI_INIT_SIMULATION_H_*/
