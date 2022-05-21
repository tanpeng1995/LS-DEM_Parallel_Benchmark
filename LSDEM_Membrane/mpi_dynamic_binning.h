/*
 * MPI_DYNAMIC_BINNING_H_
 * Created on: June 18, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_DYNAMIC_BINNING_H_
#define MPI_DYNAMIC_BINNING_H_
#include "mpi_helper_function.h"
#include "mpi_particle_migration.h"
void reassign_particles(int rank, string indir){
  if(USED > 1){
    //if(rank < USED){

      minX = INFINITY; maxX = -INFINITY;
      minY = INFINITY; maxY = -INFINITY;
      minZ = INFINITY; maxZ = -INFINITY;

      for (int i = 0; i < num_grains; ++i) {
        if(belongToThisRank[i] == 1){
          minX = min(minX, grainsWorld[i]->getPosition()(0));
          maxX = max(maxX, grainsWorld[i]->getPosition()(0));
          minY = min(minY, grainsWorld[i]->getPosition()(1));
          maxY = max(maxY, grainsWorld[i]->getPosition()(1));
          minZ = min(minZ, grainsWorld[i]->getPosition()(2));
          maxZ = max(maxZ, grainsWorld[i]->getPosition()(2));
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
      //below 3 lines is optional, if totalBinsInDimX < dimX, which means some processors are not being used.
      //if(totalBinsInDimX < dimX) totalBinsInDimX = dimX;
      //if(totalBinsInDimY < dimY) totalBinsInDimY = dimY;
      //if(totalBinsInDimZ < dimZ) totalBinsInDimZ = dimZ;
      binsInDimX = ceil((double)totalBinsInDimX/dimX) + 2; //number of bins in X dimension for a processor, including 2 ghost
      binsInDimY = ceil((double)totalBinsInDimY/dimY) + 2; //number of bins in Y dimension for a processor, including 2 ghost
      binsInDimZ = ceil((double)totalBinsInDimZ/dimZ) + 2; //number of bins in Z dimension for a processor, including 2 ghost
      totalBinsPerProc = binsInDimX*binsInDimY*binsInDimZ; //number of bins in a processor
      //cout<<"totalBinsInDimY = "<<totalBinsInDimY<<" binsInDimY = "<<binsInDimY<<" dimY = "<<dimY<<endl;
      initMPIBuff_dynamic(num_grains, SIZE, totalBinsPerProc, binsInDimX, binsInDimY, binsInDimZ); //linked list structures are re-initialzed, except for belongToThisRank
      //bins[_bin_id] is clear, bins.resize(totalBinsPerProc);
      //grainIdList[i] = -1;
      //binIdList[i] = -1;
      //firstGrainIdList[i] = -1;
      int j,_idx,_idy,_idz,_bin_id,_proc_id;
      for(int i = 0; i < num_grains; ++i){
        find_pos(grainsWorld[i]->getPosition(), _idx, _idy, _idz);
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
            if(_proc_id >= USED){
              cout<<"_id = "<<i<<" _proc_id = "<<_proc_id<<endl;
              cout<<"ERROR IN MPI_DYNAMIC_BINNING_H_, CHECK BOUNDARIES."<<endl;
              cout<<"POSSIBLE REASON: GRAINS FALL EXACTLY ON THE BOUNDARY. PRINT INFO TO DEBUG"<<endl;
              cout<<"OFTEN IN GRAIN3D_H_. FUNCTION: takeTimestep()"<<endl;
            }else if(_proc_id < 0){
              cout<<"_id = "<<i<<" _proc_id = "<<_proc_id<<endl;
              cout<<"ERROR IN MPI_DYNAMIC_BINNING_H_, CHECK SegFault."<<endl;
            }
            // if rank != _proc_id means need to migrate particle
            packer[_proc_id].emplace_back(i, grainsWorld[i]->getNodeShears().size(), grainsWorld[i]->getQuat(),
              grainsWorld[i]->getPosition(), grainsWorld[i]->getVelocity(), grainsWorld[i]->getOmega());
            packerNodeShears[_proc_id].insert(packerNodeShears[_proc_id].end(),
              grainsWorld[i]->getNodeShears().begin(),grainsWorld[i]->getNodeShears().end());
            packerNodeContact[_proc_id].insert(packerNodeContact[_proc_id].end(),
              grainsWorld[i]->getNodeContact().begin(),grainsWorld[i]->getNodeContact().end());
            packerNodeNormals[_proc_id].insert(packerNodeNormals[_proc_id].end(),
              grainsWorld[i]->getNodeNormals().begin(),grainsWorld[i]->getNodeNormals().end());

            belongToThisRank[i] = 0;
          }
        }else if(belongToThisRank[i] == 2){
          if(!ghost_bin(_idx, _idy, _idz, rank)){
            belongToThisRank[i] = 0;
          }
        }
      }
      //particle migration.
      grainsExchange(indir); // mpi_particle_migration.
    //}
  }
}
#endif /*MPI_DYNAMIC_BINNING_H_*/
