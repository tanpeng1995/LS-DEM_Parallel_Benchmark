#ifndef MPI_DYNAMIC_BINNING_HPP_
#define MPI_DYNAMIC_BINNING_HPP_
#include "mpi_helper_function.hpp"
#include "mpi_move_grain.hpp"
void reassign_particles(int rank, string indir){
  if(USED > 1){
    if(rank < USED){

      minX = INFINITY; maxX = -INFINITY;
      minY = INFINITY; maxY = -INFINITY;
      minZ = INFINITY; maxZ = -INFINITY;

      for (int i = 0; i < num_grains; ++i) {
        if(belongToThisRank[i] == 1){
          minX = min(minX, grainsWorld[i].getPosition()(0));
          maxX = max(maxX, grainsWorld[i].getPosition()(0));
          minY = min(minY, grainsWorld[i].getPosition()(1));
          maxY = max(maxY, grainsWorld[i].getPosition()(1));
          minZ = min(minZ, grainsWorld[i].getPosition()(2));
          maxZ = max(maxZ, grainsWorld[i].getPosition()(2));
        }
      }
      minX = floor(minX/CUTOFF)*CUTOFF;
      maxX = ceil(maxX/CUTOFF)*CUTOFF;
      minY = floor(minY/CUTOFF)*CUTOFF;
      maxY = ceil(maxY/CUTOFF)*CUTOFF;
      minZ = floor(minZ/CUTOFF)*CUTOFF;
      maxZ = ceil(maxZ/CUTOFF)*CUTOFF;

      MPI_Allreduce(MPI_IN_PLACE, &minX, 1, MPI_DOUBLE, MPI_MIN, CART_COMM);
      MPI_Allreduce(MPI_IN_PLACE, &maxX, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
      MPI_Allreduce(MPI_IN_PLACE, &minY, 1, MPI_DOUBLE, MPI_MIN, CART_COMM);
      MPI_Allreduce(MPI_IN_PLACE, &maxY, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
      MPI_Allreduce(MPI_IN_PLACE, &minZ, 1, MPI_DOUBLE, MPI_MIN, CART_COMM);
      MPI_Allreduce(MPI_IN_PLACE, &maxZ, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);

      worldWidth  = maxX-minX;
      worldLength = maxY-minY;
      worldHeight = maxZ-minZ;

      freeMPIBuff_dynamic();
      //during the simulation, worldWidth, worldLength and worldHeight would change
      totalBinsInDimX = ceil(worldWidth/CUTOFF); // total number of bins in X dimension
      totalBinsInDimY = ceil(worldLength/CUTOFF); // total number of bins in Y dimension
      totalBinsInDimZ = ceil(worldHeight/CUTOFF); // total number of bins in Z dimension
      binsInDimX = ceil((double)totalBinsInDimX/dimX) + 2; //number of bins in X dimension for a processor, including 2 ghost
      binsInDimY = ceil((double)totalBinsInDimY/dimY) + 2; //number of bins in Y dimension for a processor, including 2 ghost
      binsInDimZ = ceil((double)totalBinsInDimZ/dimZ) + 2; //number of bins in Z dimension for a processor, including 2 ghost
      totalBinsPerProc = binsInDimX*binsInDimY*binsInDimZ; //number of bins in a processor
      initMPIBuff_dynamic(num_grains, SIZE, totalBinsPerProc, binsInDimX, binsInDimY, binsInDimZ); //linked list structures are re-initialzed, except for belongToThisRank
      int j,_idx,_idy,_idz,_bin_id,_proc_id;
      for(int i = 0; i < num_grains; ++i){
        find_pos(grainsWorld[i].getPosition(), _idx, _idy, _idz);
        _proc_id = to_proc_id(_idx, _idy, _idz);
        _bin_id = to_bin_id(_idx, _idy, _idz);
        if(belongToThisRank[i] == 1){
          if(rank == _proc_id){
            //only need to reassign grainsWorld[i] in current rank
            bins[_bin_id].insert(i);
            binIdList[i] = _bin_id;
            if(firstGrainIdList[_bin_id] == -1){
              firstGrainIdList[_bin_id] = i;
            }else{
              j = firstGrainIdList[_bin_id];
              while(grainIdList[j] != -1) j = grainIdList[j];
              grainIdList[j] = i;
            }
          }else{
            //enum DIRECTIONS {LEFT, RIGHT, FRONT, BACK, DOWN, UP};
            //move to the left
            data_t migrated_grain = data_t(i,grainsWorld[i].getQuat(),
              grainsWorld[i].getPosition(),grainsWorld[i].getVelocity(),grainsWorld[i].getOmega());
            if(_proc_id%dimX < rank%dimX)      grain_to_be_send[LEFT].push_back(migrated_grain);
            else if(_proc_id%dimX > rank%dimX) grain_to_be_send[RIGHT].push_back(migrated_grain);
            else if(_proc_id%(dimX*dimY)/dimX < rank%(dimX*dimY)/dimX) grain_to_be_send[FRONT].push_back(migrated_grain);
            else if(_proc_id%(dimX*dimY)/dimX > rank%(dimX*dimY)/dimX) grain_to_be_send[BACK].push_back(migrated_grain);
            else if(_proc_id/(dimX*dimY) < rank/(dimX*dimY)) grain_to_be_send[DOWN].push_back(migrated_grain);
            else if(_proc_id/(dimX*dimY) > rank/(dimX*dimY)) grain_to_be_send[UP].push_back(migrated_grain);

            belongToThisRank[i] = 0;
          }
        }else if(belongToThisRank[i] == 2){
          if(!ghost_bin(_idx, _idy, _idz, rank)){
            belongToThisRank[i] = 0;
          }
        }
      }
      //particle migration.
      grainsExchange(rank, indir); // mpi_particle_migration.
    }
  }
}
#endif /*MPI_DYNAMIC_BINNING_HPP_*/
