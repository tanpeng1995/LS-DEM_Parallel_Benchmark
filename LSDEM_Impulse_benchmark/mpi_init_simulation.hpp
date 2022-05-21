#ifndef MPI_INIT_SIMULATION_H_
#define MPI_INIT_SIMULATION_H_
#include "mpi_helper_function.hpp"
#include "Grain3d.hpp"
#include "WallPlane.hpp"

void init_simulation(int rank, int num_procs, string indir){
  //determine used processor and build new communicator group.
  //num_procs-1 because need to leave one processor doint boundary staff
  determine_proc(num_grains, num_procs, USED, dimX, dimY, dimZ);
  sort(dimX, dimY, dimZ); //so that dimZ >= dimX >= dimY
  if(rank == 0) cout<<"dimZ: "<<dimZ<<" dimY: "<<dimY<<" dimX: "<<dimX<<endl;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);

  int ranks[USED];
  for(int i = 0; i < USED; ++i) ranks[i] = i;
  MPI_Group_incl(world_group, USED, ranks, &new_group);
  MPI_Comm_create(MPI_COMM_WORLD, new_group, &USING_COMM); // create communicator USING_COMM this is we actually use for case without boundary.
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

  const int items2 = 3;
  int lengths2[3] = {1, 4, 3};
  MPI_Datatype types2[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offsets2[3];
  offsets2[0] = offsetof(output_data, _id);
  offsets2[1] = offsetof(output_data, _quat);
  offsets2[2] = offsetof(output_data, _position);
  MPI_Type_create_struct(items2, lengths2, offsets2, types2, &OUTPUT_DATA);
  MPI_Type_commit(&OUTPUT_DATA);
  MPI_Type_contiguous(SIZE, DATA, &BIN);
  MPI_Type_commit(&BIN);

  const int items3 = 4;
  int lengths3[4] = {1, 1, 3, 3};
  MPI_Datatype types3[4] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offsets3[4];
  offsets3[0] = offsetof(vel_data, _id);
  offsets3[1] = offsetof(vel_data, _rank);
  offsets3[2] = offsetof(vel_data, _velocity);
  offsets3[3] = offsetof(vel_data, _omegaGlobal);
  MPI_Type_create_struct(items3, lengths3, offsets3, types3, &VELOCITY);
  MPI_Type_commit(&VELOCITY);

  const int items4 = 8;
  int lengths4[8] = {1, 1, 3, 3, 3, 3, 1, 1};
  MPI_Datatype types4[8] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offsets4[8];
  offsets4[0] = offsetof(ContactInfo, _master);
  offsets4[1] = offsetof(ContactInfo, _slave);
  offsets4[2] = offsetof(ContactInfo, _masterR);
  offsets4[3] = offsetof(ContactInfo, _slaveR);
  offsets4[4] = offsetof(ContactInfo, _normal);
  offsets4[5] = offsetof(ContactInfo, _velocity);
  offsets4[6] = offsetof(ContactInfo, _normalV);
  offsets4[7] = offsetof(ContactInfo, _energy);
  MPI_Type_create_struct(items4, lengths4, offsets4, types4, &CONTACT);
  MPI_Type_commit(&CONTACT);

  if(rank < USED){
    vector<Vector3d> inputPos;
    vector<Vector4d> inputRot;

    inputPos = readPositionFile(indir + "InitState/positions_"+testName+".dat", target_grains);
    inputRot = readQuaternionFile(indir + "InitState/rotations_"+testName+".dat", target_grains);

    grainsWorld = std::make_unique<Grain3d[]>(num_grains);
    grainIdList = std::make_unique<int[]>(num_grains);
    binIdList = std::make_unique<int[]>(num_grains);
    belongToThisRank = std::make_unique<int[]>(num_grains);
    if(rank == 0) cout<<"finish constructing linked-list..."<<endl;

    wallPlanes = std::make_unique<WallPlane[]>(6);
    wallPlanes[LEFT]  = WallPlane(Vector3d(1.,0.,0.), Vector3d(0.,0.,0.), wallStiffness);
    wallPlanes[RIGHT] = WallPlane(Vector3d(-1.,0.,0.), Vector3d(worldWidth,0.,0.), wallStiffness);
    wallPlanes[FRONT] = WallPlane(Vector3d(0.,-1.,0.), Vector3d(0.,worldLength,0.), wallStiffness);
    wallPlanes[BACK]  = WallPlane(Vector3d(0.,1.,0.), Vector3d(0.,0.,0.), wallStiffness);
    wallPlanes[DOWN]  = WallPlane(Vector3d(0.,0.,1.), Vector3d(0.,0.,0.), wallStiffness);
    wallPlanes[UP]    = WallPlane(Vector3d(0.,0.,-1.), Vector3d(0.,0.,worldHeight), wallStiffness);
    if(rank == 0) cout<<"finish buiding wallPlane and wallCap..."<<endl;

    maxX = worldWidth; maxY = worldLength; maxZ = worldHeight;
    minX = 0.; minY = 0.; minZ = 0.;

    worldWidth  = maxX-minX;
    worldLength = maxY-minY;
    worldHeight = maxZ-minZ;

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

    minMass = std::numeric_limits<double>::max();
    maxMass = std::numeric_limits<double>::min();
    int _idx, _idy, _idz, _bin_id, _proc_id;
    int Nx = worldWidth / targetWidth;
    int Ny = worldLength / targetLength;
    int Nz = worldHeight / targetHeight;
    for(int iz = 0; iz < Nz; ++iz){
      for(int iy = 0; iy < Ny; ++iy){
        for(int ix = 0; ix < Nx; ++ix){
          Vector3d disp(ix*targetWidth,iy*targetLength,iz*targetHeight);
          int disp_i = (Nx*Ny*iz + Nx*iy + ix)*target_grains;
          /* BEGAIN DUPLICATING */
          for(int i = 0; i < target_grains; ++i){
            int gid = i + disp_i;
            if(rank == 0 && gid % 2000 == 0) cout<<"finish loading "<<gid<<" grains"<<endl;
            find_pos(inputPos[i]+disp, _idx, _idy, _idz);
            _proc_id = to_proc_id(_idx, _idy, _idz);
            _bin_id = to_bin_id(_idx, _idy, _idz);
            stringstream morphfile;
            morphfile << indir << "Morphologies/morph_" << i+1 << ".dat";
            if(_proc_id == rank){
              grainsWorld[gid] = generateGrainFromFile(inputPos[i]+disp, inputRot[i], morphfile.str(), grainDensity, gid);
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
              double grainMass = grainsWorld[gid].getMass();
              minMass = minMass < grainMass ? minMass : grainMass;
              maxMass = maxMass > grainMass ? maxMass : grainMass;
            }else if(ghost_bin(_idx, _idy, _idz, rank)){
              grainsWorld[gid] = generateGrainFromFile(inputPos[i]+disp, inputRot[i], morphfile.str(), grainDensity, gid);
              grainsWorld[gid].setRealGrain(true);
              belongToThisRank[gid] = 2;
            }else{
              grainsWorld[gid] = generateSimpleGrain3d(inputPos[i]+disp, inputRot[i], morphfile.str(), grainDensity, gid);
              grainsWorld[gid].setRealGrain(false);
              belongToThisRank[gid] = 0;
            }
          }
          /* END DUPLICATING */
        }
      }
    }

    /* broadcast minMass and maxMass */
    MPI_Allreduce(MPI_IN_PLACE, &minMass, 1, MPI_DOUBLE, MPI_MIN, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &maxMass, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    if(massScaling){
      for(int i = 0; i < num_grains; ++i){
        double grainMass = grainsWorld[i].getMass();
        if(grainMass < maxMass / 10.){
          grainMass = maxMass / 10.;
          grainsWorld[i].changeMass(grainMass);
        }
      }
    }

  }
}
#endif /*MPI_INIT_SIMULATION_HPP_*/
