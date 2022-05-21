
/*
 * MPI_INIT_SIMULATION_H_
 *
 * Created on: June 20, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_INIT_SIMULATION_H_
#define MPI_INIT_SIMULATION_H_
#include "mpi_helper_function.h"
#include "mpi_grains.h"

void init_simulation(int rank, int num_procs, string indir, FILE * wallfile, FILE * wallfile2){
  //determine used processor and build new communicator group.
  //num_procs-1 because need to leave one processor doint boundary staff
  determine_proc(num_grains, num_procs, USED, dimX, dimY, dimZ);
  sort(dimX, dimY, dimZ); //so that dimZ >= dimX >= dimY
  if(rank == 0) cout<<"dimZ: "<<dimZ<<" dimY: "<<dimY<<" dimX: "<<dimX<<endl;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);

  int ranks[USED];
  for(int i = 0; i < USED; ++i) ranks[i] = i;
  MPI_Group_incl(world_group, USED, ranks, &new_group);
  MPI_Comm_create(MPI_COMM_WORLD, new_group, &USING_COMM);
  if(rank == 0) cout<<"finish building USING_COMM..."<<endl;

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

  const int items4 = 3;
  int lengths4[3] = {1, 1, 3};
  MPI_Datatype types4[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};
  MPI_Aint offsets4[3];
  offsets4[0] = offsetof(MPIBalls, _rank);
  offsets4[1] = offsetof(MPIBalls, _id);
  offsets4[2] = offsetof(MPIBalls, _force);
  MPI_Type_create_struct(items4, lengths4, offsets4, types4, &MPI_BALL);
  MPI_Type_commit(&MPI_BALL);
  if(rank == 0) cout<<"finish creating MPI datatype..."<<endl;

  MPI_Type_contiguous(SIZE, DATA, &BIN);
  MPI_Type_commit(&BIN);

  if(rank < USED){
    vector<Vector3d> inputPos;
    vector<Vector4d> inputRot;

    inputPos = readPositionFile(indir + InitState + "/positions_"+testName+".dat", num_grains);
    inputRot = readQuaternionFile(indir + InitState + "/rotations_"+testName+".dat", num_grains);

    grainsWorld.resize(num_grains);
    grainIdList = std::make_unique<int[]>(num_grains);
    binIdList = std::make_unique<int[]>(num_grains);
    belongToThisRank = std::make_unique<int[]>(num_grains);
    if(rank == 0) cout<<"finish constructing linked-list..."<<endl;

    wallPlane = std::make_unique<WallPlane> (Vector3d(0.,0.,1.), Vector3d(0.,0.,0.), kn, wallFriction, initRadius, num_grains);
    wallCap = std::make_unique<WallCap>(Vector3d(0.,0.,-1.), Vector3d(0.,0.,initHeight), kn, wallFriction, num_grains+1, gDamping, capHeight, capRadius, capDensity);
    //wallCylinder = std::make_unique<WallCylinder>(initHeight, initRadius, kn, ks, mu, num_grains+2);
    wallBalls = std::make_unique<WallBalls> (initHeight, initRadius, ballRadius, 1.0*pressure, wallDensity, wallToGrainStiffness, kmn, kms, gDamping);
    externalWall = std::make_unique<WallBalls>(initHeight, initRadius+2.0*ballRadius, 1.1*ballRadius, 0.2*pressure, wallDensity, wallToGrainStiffness, kmn, kms, gDamping);
    if(rank == 0) cout<<"finish buiding wallPlane and wallCap..."<<endl;

    //find domain boundary
    for(int i = 0; i < num_grains; ++i){
      minX    = min(minX, inputPos[i](0));
      maxX    = max(maxX, inputPos[i](0));
      minY    = min(minY, inputPos[i](1));
      maxY    = max(maxY, inputPos[i](1));
      minZ    = min(minZ, inputPos[i](2));
      maxZ    = max(maxZ, inputPos[i](2));
    }
    //pad minX, maxX, minY, maxY, minZ, maxZ, make min smaller and max larger.
    //for every sort of boundary, this section should be modified.
    minX = floor(minX/CUTOFF)*CUTOFF;
    maxX = ceil(maxX/CUTOFF)*CUTOFF;
    minY = floor(minY/CUTOFF)*CUTOFF;
    maxY = ceil(maxY/CUTOFF)*CUTOFF;
    minZ = floor(minZ/CUTOFF)*CUTOFF;
    maxZ = ceil(maxZ/CUTOFF)*CUTOFF;

    worldWidth = maxX-minX;
    worldLength = maxY-minY;
    worldHeight = maxZ-minZ;
    stressVoigt<<0., 0., 0., 0., 0., 0.;

    // if rank is normal processor
    int dims[3] = {dimZ, dimY, dimX};
    MPI_Dims_create(USED, 3, dims);
    int periods[3] = {false, false, false};
    int reorder = false;
    MPI_Cart_create(USING_COMM, 3, dims, periods, reorder, &CART_COMM);
    //if this line put outside the loop, will report error, because USING_COMM only includes rank < USED.
    MPI_Cart_shift(CART_COMM, 0, 1, &neighbor[DOWN], &neighbor[UP]);
    MPI_Cart_shift(CART_COMM, 1, 1, &neighbor[FRONT], &neighbor[BACK]);
    MPI_Cart_shift(CART_COMM, 2, 1, &neighbor[LEFT], &neighbor[RIGHT]);

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
    if(rank == 0) cout<<"finish creating bin, blocks and allocating memories..."<<endl;

    //assign grains
    int _idx, _idy, _idz, _bin_id, _proc_id;
    for(int i = 0; i < num_grains; ++i){
      if(rank == 0 && i % 1000 == 0) cout<<"finish loading "<<i<<" grains"<<endl;
      find_pos(inputPos[i], _idx, _idy, _idz);
      _proc_id = to_proc_id(_idx, _idy, _idz);
      _bin_id = to_bin_id(_idx, _idy, _idz);
      stringstream morphfile;
      morphfile << indir << "Morphologies/morph_" << i+1 << ".dat";
      if(_proc_id == rank){
        grainsWorld[i] = generateGrainFromFile(inputPos[i], inputRot[i], morphfile.str(), i);
        grainsWorld[i]->changeDensity(grainDensity);
        grainsWorld[i]->setRealGrain(true);
        bins[_bin_id].insert(i);
        binIdList[i] = _bin_id;

        if(firstGrainIdList[_bin_id] == -1){
          firstGrainIdList[_bin_id] = i;
        }else{
          int j = firstGrainIdList[_bin_id];
          while(grainIdList[j] != -1) j = grainIdList[j];
          grainIdList[j] = i;
        }
        belongToThisRank[i] = 1;
        minMass = min(minMass, grainsWorld[i]->getMass());
        maxMass = max(maxMass, grainsWorld[i]->getMass());
        maxR    = max(maxR, grainsWorld[i]->getRadius());
        avgRadius += grainsWorld[i]->getRadius();
        totalMass += grainsWorld[i]->getMass();
      }else if(ghost_bin(_idx, _idy, _idz, rank)){
        grainsWorld[i] = generateGrainFromFile(inputPos[i], inputRot[i], morphfile.str(), i);
        grainsWorld[i]->changeDensity(grainDensity);
        grainsWorld[i]->setRealGrain(true);
        belongToThisRank[i] = 2;
      }else{
        grainsWorld[i] = std::make_unique<Grain3d>(inputPos[i], inputRot[i], i);
        grainsWorld[i]->setRealGrain(false);
        belongToThisRank[i] = 0;
      }
    }

    fprintf(wallfile, "%d\n", wallBalls->getNumWalls());
    fprintf(wallfile, "%.4f\n", wallBalls->getBallRadius());
    fprintf(wallfile2, "%d\n", externalWall->getNumWalls());
    fprintf(wallfile2, "%.4f\n", externalWall->getBallRadius());
    fflush(wallfile);
    fflush(wallfile2);

    MPI_Allreduce(MPI_IN_PLACE, &minMass, 1, MPI_DOUBLE, MPI_MIN, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &maxMass, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &maxR, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &avgRadius, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &totalMass, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
    avgRadius = avgRadius/num_grains;

    if(massScaling){ minMass = maxMass / massScalingCoefficient; }
    if(massScaling){
      for(int i = 0; i < num_grains; ++i){
        if(belongToThisRank[i] != 0){
          if(grainsWorld[i]->getMass() < minMass){
            grainsWorld[i]->multiplyMomentInertia(minMass/grainsWorld[i]->getMass());
            grainsWorld[i]->changeMass(minMass);
          }
        }
      }
    }
    /*end massScaling*/
  }
}
#endif /*MPI_INIT_SIMULATION_H_*/
