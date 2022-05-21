/*
 * MPI_PARTICLE_MIGRATION_H_
 * Created on: June 13, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_PARTICLE_MIGRATION_H_
#define MPI_PARTICLE_MIGRATION_H_
#include "mpi_helper_function.h"
#include "mpi_grains.h"
void grainsExchange(string indir){
  for(int i = 0; i < USED; ++i){
    countSendGrain[i] = 0;
    dispSendGrain[i]  = 0;
    countRecvGrain[i] = 0;
    dispRecvGrain[i]  = 0;

    countSendVector3d[i] = 0;
    dispSendVector3d[i]  = 0;
    countRecvVector3d[i] = 0;
    dispRecvVector3d[i]  = 0;

    countSendInt[i] = 0;
    dispSendInt[i]  = 0;
    countRecvInt[i] = 0;
    dispRecvInt[i]  = 0;
  }
  totalGrainRecv = 0;
  totalNodeRecv  = 0;

  allToAllSendGrains.clear();
  allToAllSendGrains.shrink_to_fit();
  allTOAllRecvGrains.clear();
  allTOAllRecvGrains.shrink_to_fit();
  allToAllSendNodeShears.clear();
  allToAllSendNodeShears.shrink_to_fit();
  allToAllRecvNodeShears.clear();
  allToAllRecvNodeShears.shrink_to_fit();
  allToAllSendNodeContact.clear();
  allToAllSendNodeContact.shrink_to_fit();
  allToAllRecvNodeContact.clear();
  allToAllRecvNodeContact.shrink_to_fit();
  allToAllSendNodeNormals.clear();
  allToAllSendNodeNormals.shrink_to_fit();
  allToAllRecvNodeNormals.clear();
  allToAllRecvNodeNormals.shrink_to_fit();

  //build all send buffer
  for(int i = 0; i < USED; ++i){
    countSendGrain[i] = packer[i].size();
    countSendVector3d[i] = packerNodeShears[i].size()*3;
    countSendInt[i] = packerNodeShears[i].size();

    allToAllSendGrains.insert(allToAllSendGrains.end(),packer[i].begin(),packer[i].end());
    allToAllSendNodeShears.insert(allToAllSendNodeShears.end(),packerNodeShears[i].begin(),packerNodeShears[i].end());
    allToAllSendNodeContact.insert(allToAllSendNodeContact.end(),packerNodeContact[i].begin(),packerNodeContact[i].end());
    allToAllSendNodeNormals.insert(allToAllSendNodeNormals.end(),packerNodeNormals[i].begin(),packerNodeNormals[i].end());

    packer[i].clear();
    packer[i].shrink_to_fit();
    packerNodeShears[i].clear();
    packerNodeShears[i].shrink_to_fit();
    packerNodeContact[i].clear();
    packerNodeContact[i].shrink_to_fit();
    packerNodeNormals[i].clear();
    packerNodeNormals[i].shrink_to_fit();

    if(i == 0){
      dispSendGrain[i] = 0;
      dispSendVector3d[i] = 0;
      dispSendInt[i] = 0;
    }else{
      dispSendGrain[i] = countSendGrain[i-1]+dispSendGrain[i-1];
      dispSendVector3d[i] = countSendVector3d[i-1]+dispSendVector3d[i-1];
      dispSendInt[i]   = countSendInt[i-1]+dispSendInt[i-1];
    }
  }
  //alltoall
  MPI_Alltoall(countSendGrain,1,MPI_INT,countRecvGrain,1,MPI_INT,CART_COMM);
  MPI_Alltoall(countSendVector3d,1,MPI_INT,countRecvVector3d,1,MPI_INT,CART_COMM);
  MPI_Alltoall(countSendInt,1,MPI_INT,countRecvInt,1,MPI_INT,CART_COMM);
  //build all recv buffer
  for(int i = 0; i < USED; ++i){
    totalGrainRecv += countRecvGrain[i];
    totalNodeRecv += countRecvInt[i];

    if(i == 0){
      dispRecvGrain[i] = 0;
      dispRecvVector3d[i] = 0;
      dispRecvInt[i] = 0;
    }else{
      dispRecvGrain[i] = countRecvGrain[i-1]+dispRecvGrain[i-1];
      dispRecvVector3d[i] = countRecvVector3d[i-1]+dispRecvVector3d[i-1];
      dispRecvInt[i] = countRecvInt[i-1]+dispRecvInt[i-1];
    }
  }

  allTOAllRecvGrains.resize(totalGrainRecv);
  allToAllRecvNodeShears.resize(totalNodeRecv);
  allToAllRecvNodeContact.resize(totalNodeRecv);
  allToAllRecvNodeNormals.resize(totalNodeRecv);
  //MPI_Alltoallv(buffer_send, counts_send, displacements_send, MPI_Datatype, buffer_recv, counts_recv, displacements_recv, MPI_Datatype, MPI_COMM_WORLD);
  // counts_send: # of items sent to each rank;
  // displacements_send: displacements of each items (in items not in bytes)
  MPI_Alltoallv(allToAllSendGrains.data(),countSendGrain,dispSendGrain,COMPLETE_DATA,
                allTOAllRecvGrains.data(),countRecvGrain,dispRecvGrain,COMPLETE_DATA,CART_COMM);
  MPI_Alltoallv(allToAllSendNodeShears.data(),countSendVector3d,dispSendVector3d,MPI_DOUBLE,
                allToAllRecvNodeShears.data(),countRecvVector3d,dispRecvVector3d,MPI_DOUBLE,CART_COMM);
  MPI_Alltoallv(allToAllSendNodeContact.data(),countSendInt,dispSendInt,MPI_INT,
                allToAllRecvNodeContact.data(),countRecvInt,dispRecvInt,MPI_INT,CART_COMM);
  MPI_Alltoallv(allToAllSendNodeNormals.data(),countSendVector3d,dispSendVector3d,MPI_DOUBLE,
                allToAllRecvNodeNormals.data(),countRecvVector3d,dispRecvVector3d,MPI_DOUBLE,CART_COMM);

  vector<Vector3d> recvNodeShears;
  vector<int> recvNodeContact;
  vector<Vector3d> recvNodeNormals;
  int curr = 0, npoints = 0;
  for(int i = 0; i < totalGrainRecv; ++i){
    int gid = allTOAllRecvGrains[i]._id;
    int _idx, _idy, _idz;
    find_pos(allTOAllRecvGrains[i]._position, _idx, _idy, _idz);
    int _bin_id = to_bin_id(_idx,_idy,_idz);
    bins[_bin_id].insert(gid);

    //grainsWorld[gid].changePosition(allTOAllRecvGrains[i]._position); // update _position
    //grainsWorld[gid].changeRotation(allTOAllRecvGrains[i]._quat); // update _quat, update _rotMatrix, update _pointList
    //grainsWorld[gid].changeVelocity(allTOAllRecvGrains[i]._velocity); // update _velocity
    //grainsWorld[gid].changeOmega(allTOAllRecvGrains[i]._omega); // update _omega, update _omegaGlobal
    Vector3d position = allTOAllRecvGrains[i]._position;
    Vector3d velocity = allTOAllRecvGrains[i]._velocity;
    Vector3d omega = allTOAllRecvGrains[i]._omega;
    Vector4d quat(allTOAllRecvGrains[i]._quat[0], allTOAllRecvGrains[i]._quat[1],
                  allTOAllRecvGrains[i]._quat[2], allTOAllRecvGrains[i]._quat[3]);
    if(grainsWorld[gid].getRealGrain()){
/*
      if(belongToThisRank[gid] == 1){
        cout<<"something must be wrong with grain: "<<gid<<endl;
      }else if(belongToThisRank[gid] == 2){
        cout<<"ghost grain moves in."<<endl;
      }
*/
      grainsWorld[gid].changePosition(allTOAllRecvGrains[i]._position); // update _position
      grainsWorld[gid].changeRotation(allTOAllRecvGrains[i]._quat); // update _quat, update _rotMatrix, update _pointList
      grainsWorld[gid].changeVelocity(allTOAllRecvGrains[i]._velocity); // update _velocity
      grainsWorld[gid].changeOmega(allTOAllRecvGrains[i]._omega); // update _omega, update _omegaGlobal
    }else{
      //cout<<gid<<" will be read from file..."<<endl;
      stringstream morphfile;
      morphfile << indir << "Morphologies/morph_" << gid+1 << ".dat";
      grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), gid);
      grainsWorld[gid].changeDensity(grainDensity);
      grainsWorld[gid].setRealGrain(true);
    }

    npoints = allTOAllRecvGrains[i]._npoints;
    recvNodeShears.assign(allToAllRecvNodeShears.begin()+curr, allToAllRecvNodeShears.begin()+curr+npoints);
    recvNodeContact.assign(allToAllRecvNodeContact.begin()+curr, allToAllRecvNodeContact.begin()+curr+npoints);
    recvNodeNormals.assign(allToAllRecvNodeNormals.begin()+curr, allToAllRecvNodeNormals.begin()+curr+npoints);
    grainsWorld[gid].changeShearHist(recvNodeShears,recvNodeContact,recvNodeNormals);
    recvNodeShears.shrink_to_fit(); recvNodeShears.shrink_to_fit(); recvNodeNormals.shrink_to_fit();
    curr += npoints;

    //ADD
    if(firstGrainIdList[_bin_id] == -1){
      firstGrainIdList[_bin_id] = gid;
    }else{
      int j = firstGrainIdList[_bin_id];
      while(grainIdList[j] != -1) j = grainIdList[j];
      grainIdList[j] = gid;
    }
    binIdList[gid] = _bin_id;
    belongToThisRank[gid] = 1;
  }
  //END
}
#endif /*MPI_PARTICLE_MIGRATION_H_*/
