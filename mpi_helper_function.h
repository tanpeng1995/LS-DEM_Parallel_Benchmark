/*
 * MPI_HELPER_FUNCTION_H_
 * for simple numerical experiment on the effects binning algorithm
 * Created on: May 6, 2020
 *      Author: Peng TAN (Berkeley)
 *
 */
#ifndef MPI_HELPER_FUNCTION_H_
#define MPI_HELPER_FUNCTION_H_
#include "readParameters.h"
inline int to_bin_id(int idx, int idy, int idz){
  return (idz%(binsInDimZ-2)+1)*binsInDimX*binsInDimY + (idy%(binsInDimY-2)+1)*binsInDimX + (idx%(binsInDimX-2)+1);
}
inline int to_proc_id(int idx, int idy, int idz){
  return (idz/(binsInDimZ-2)*dimX*dimY + idy/(binsInDimY-2)*dimX + idx/(binsInDimX-2));
}
inline void find_pos(const Vector3d & position, int & idx, int & idy, int & idz){
  idx = int((position(0)-minX)/CUTOFF);
  idy = int((position(1)-minY)/CUTOFF);
  idz = int((position(2)-minZ)/CUTOFF);

  if(idx < 0 || idx >= totalBinsInDimX){ idx = idx < 0? 0 : totalBinsInDimX-1; }
  if(idy < 0 || idy >= totalBinsInDimY){ idy = idy < 0? 0 : totalBinsInDimY-1; }
  if(idz < 0 || idz >= totalBinsInDimZ){ idz = idz < 0? 0 : totalBinsInDimZ-1; }
}

inline void find_pos(const int & binNum, int & idx, int & idy, int & idz){
  idx = binNum%(binsInDimX*binsInDimY)%binsInDimX;
  idy = binNum%(binsInDimX*binsInDimY)/binsInDimX;
  idz = binNum/(binsInDimX*binsInDimY);
}

void initMPIBuff_static(const int & num_grains, const int & USED){
  for(int i = 0; i < num_grains; ++i){
    grainIdList[i] = -1;
    binIdList[i] = -1;
    belongToThisRank[i] = 0;
  }

  packer.resize(USED);
  packerNodeShears.resize(USED);
  packerNodeContact.resize(USED);
  packerNodeNormals.resize(USED);

  countSendGrain = new int[USED];
  dispSendGrain  = new int[USED];
  countRecvGrain = new int[USED];
  dispRecvGrain  = new int[USED];
  countSendVector3d = new int[USED];
  dispSendVector3d  = new int[USED];
  countRecvVector3d = new int[USED];
  dispRecvVector3d  = new int[USED];
  countSendInt = new int[USED];
  dispSendInt  = new int[USED];
  countRecvInt = new int[USED];
  dispRecvInt  = new int[USED];

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
}

void initMPIBuff_dynamic(const int & num_grains, const int & SIZE, const int & totalBinsPerProc, const int & binsInDimX, const int & binsInDimY, const int & binsInDimZ){
  for(int i = 0; i < num_grains; ++i){
    grainIdList[i] = -1;
    binIdList[i] = -1;
  }

  firstGrainIdList = std::make_unique<int[]>(totalBinsPerProc);
  bins.resize(totalBinsPerProc);

  for(int i = 0; i < totalBinsPerProc; ++i){
    firstGrainIdList[i] = -1;
    bins[i].clear();
  }

  MPI_Alloc_mem(SIZE*(binsInDimY-2)*(binsInDimZ-2)*sizeof(data_t),MPI_INFO_NULL,&sendToRight);
  for(int i = 0; i < (binsInDimY-2)*(binsInDimZ-2); ++i){
    for(int j = 0; j < SIZE; ++j){
      sendToRight[i][j] = data_t();
    }
  }
  MPI_Alloc_mem(SIZE*(binsInDimY-2)*(binsInDimZ-2)*sizeof(data_t),MPI_INFO_NULL,&sendToLeft);
  for(int i = 0; i < (binsInDimY-2)*(binsInDimZ-2); ++i){
    for(int j = 0; j < SIZE; ++j){
      sendToLeft[i][j] = data_t();
    }
  }
  MPI_Alloc_mem(SIZE*binsInDimX*(binsInDimZ-2)*sizeof(data_t),MPI_INFO_NULL,&sendToBack);
  for(int i = 0; i < binsInDimX*(binsInDimZ-2); ++i){
    for(int j = 0; j < SIZE; ++j){
      sendToBack[i][j] = data_t();
    }
  }
  MPI_Alloc_mem(SIZE*binsInDimX*(binsInDimZ-2)*sizeof(data_t),MPI_INFO_NULL,&sendToFront);
  for(int i = 0; i < binsInDimX*(binsInDimZ-2); ++i){
    for(int j = 0; j < SIZE; ++j){
      sendToFront[i][j] = data_t();
    }
  }
  MPI_Alloc_mem(SIZE*binsInDimX*binsInDimY*sizeof(data_t),MPI_INFO_NULL,&sendToDown);
  for(int i = 0; i < binsInDimX*binsInDimY; ++i){
    for(int j = 0; j < SIZE; ++j){
      sendToDown[i][j] = data_t();
    }
  }
  MPI_Alloc_mem(SIZE*binsInDimX*binsInDimY*sizeof(data_t),MPI_INFO_NULL,&sendToUp);
  for(int i = 0; i < binsInDimX*binsInDimY; ++i){
    for(int j = 0; j < SIZE; ++j){
      sendToUp[i][j] = data_t();
    }
  }
  MPI_Alloc_mem((binsInDimY-2)*(binsInDimZ-2)*sizeof(int),MPI_INFO_NULL,&grainNumSendToRight);
  for(int i = 0; i < (binsInDimY-2)*(binsInDimZ-2); ++i) grainNumSendToRight[i] = 0;

  MPI_Alloc_mem((binsInDimY-2)*(binsInDimZ-2)*sizeof(int),MPI_INFO_NULL,&grainNumSendToLeft);
  for(int i = 0; i < (binsInDimY-2)*(binsInDimZ-2); ++i) grainNumSendToLeft[i] = 0;

  MPI_Alloc_mem(binsInDimX*(binsInDimZ-2)*sizeof(int),MPI_INFO_NULL,&grainNumSendToBack);
  for(int i = 0; i < binsInDimX*(binsInDimZ-2); ++i) grainNumSendToBack[i] = 0;

  MPI_Alloc_mem(binsInDimX*(binsInDimZ-2)*sizeof(int),MPI_INFO_NULL,&grainNumSendToFront);
  for(int i = 0; i < binsInDimX*(binsInDimZ-2); ++i) grainNumSendToFront[i] = 0;

  MPI_Alloc_mem(binsInDimX*binsInDimY*sizeof(int),MPI_INFO_NULL,&grainNumSendToDown);
  for(int i = 0; i < binsInDimX*binsInDimY; ++i) grainNumSendToDown[i] = 0;

  MPI_Alloc_mem(binsInDimX*binsInDimY*sizeof(int),MPI_INFO_NULL,&grainNumSendToUp);
  for(int i = 0; i < binsInDimX*binsInDimY; ++i) grainNumSendToDown[i] = 0;

  MPI_Alloc_mem(SIZE*(binsInDimY-2)*(binsInDimZ-2)*sizeof(data_t),MPI_INFO_NULL,&recvFromRight);
  for(int i = 0; i < (binsInDimY-2)*(binsInDimZ-2); ++i){
    for(int j = 0; j < SIZE; ++j){
      recvFromRight[i][j] = data_t();
    }
  }
  MPI_Alloc_mem(SIZE*(binsInDimY-2)*(binsInDimZ-2)*sizeof(data_t),MPI_INFO_NULL,&recvFromLeft);
  for(int i = 0; i < (binsInDimY-2)*(binsInDimZ-2); ++i){
    for(int j = 0; j < SIZE; ++j){
      recvFromLeft[i][j] = data_t();
    }
  }
  MPI_Alloc_mem(SIZE*binsInDimX*(binsInDimZ-2)*sizeof(data_t),MPI_INFO_NULL,&recvFromBack);
  for(int i = 0; i < binsInDimX*(binsInDimZ-2); ++i){
    for(int j = 0; j < SIZE; ++j){
      recvFromBack[i][j] = data_t();
    }
  }
  MPI_Alloc_mem(SIZE*binsInDimX*(binsInDimZ-2)*sizeof(data_t),MPI_INFO_NULL,&recvFromFront);
  for(int i = 0; i < binsInDimX*(binsInDimZ-2); ++i){
    for(int j = 0; j < SIZE; ++j){
      recvFromFront[i][j] = data_t();
    }
  }
  MPI_Alloc_mem(SIZE*binsInDimX*binsInDimY*sizeof(data_t),MPI_INFO_NULL,&recvFromDown);
  for(int i = 0; i < binsInDimX*binsInDimY; ++i){
    for(int j = 0; j < SIZE; ++j){
      recvFromDown[i][j] = data_t();
    }
  }
  MPI_Alloc_mem(SIZE*binsInDimX*binsInDimY*sizeof(data_t),MPI_INFO_NULL,&recvFromUp);
  for(int i = 0; i < binsInDimX*binsInDimY; ++i){
    for(int j = 0; j < SIZE; ++j){
      recvFromUp[i][j] = data_t();
    }
  }
  MPI_Alloc_mem((binsInDimY-2)*(binsInDimZ-2)*sizeof(int),MPI_INFO_NULL,&grainNumRecvFromRight);
  for(int i = 0; i < (binsInDimY-2)*(binsInDimZ-2); ++i) grainNumRecvFromRight[i] = 0;

  MPI_Alloc_mem((binsInDimY-2)*(binsInDimZ-2)*sizeof(int),MPI_INFO_NULL,&grainNumRecvFromLeft);
  for(int i = 0; i < (binsInDimY-2)*(binsInDimZ-2); ++i) grainNumRecvFromLeft[i] = 0;

  MPI_Alloc_mem(binsInDimX*(binsInDimZ-2)*sizeof(int),MPI_INFO_NULL,&grainNumRecvFromBack);
  for(int i = 0; i < binsInDimX*(binsInDimZ-2); ++i) grainNumRecvFromBack[i] = 0;

  MPI_Alloc_mem(binsInDimX*(binsInDimZ-2)*sizeof(int),MPI_INFO_NULL,&grainNumRecvFromFront);
  for(int i = 0; i < binsInDimX*(binsInDimZ-2); ++i) grainNumRecvFromFront[i] = 0;

  MPI_Alloc_mem(binsInDimX*binsInDimY*sizeof(int),MPI_INFO_NULL,&grainNumRecvFromDown);
  for(int i = 0; i < binsInDimX*binsInDimY; ++i) grainNumRecvFromDown[i] = 0;

  MPI_Alloc_mem(binsInDimX*binsInDimY*sizeof(int),MPI_INFO_NULL,&grainNumRecvFromUp);
  for(int i = 0; i < binsInDimX*binsInDimY; ++i) grainNumRecvFromUp[i] = 0;
}

void freeMPIBuff_static(){
  MPI_Comm_free(&USING_COMM);
  //MPI_Comm_free(&BOUNDARY_COMM);
  MPI_Comm_free(&CART_COMM);
  MPI_Type_free(&DATA);
  MPI_Type_free(&COMPLETE_DATA);
  MPI_Type_free(&OUTPUT_DATA);
  MPI_Type_free(&BIN);
  MPI_Group_free(&world_group);
  MPI_Group_free(&new_group);
  //MPI_Group_free(&boundary_group);
  bins.clear();
  bins.shrink_to_fit();
  packer.clear();
  packer.shrink_to_fit();
  packerNodeShears.clear();
  packerNodeShears.shrink_to_fit();
  packerNodeContact.clear();
  packerNodeContact.shrink_to_fit();
  packerNodeNormals.clear();
  packerNodeNormals.shrink_to_fit();
  delete [] countSendGrain;
  delete [] dispSendGrain;
  delete [] countRecvGrain;
  delete [] dispRecvGrain;
  delete [] countSendVector3d;
  delete [] dispSendVector3d;
  delete [] countRecvVector3d;
  delete [] dispRecvVector3d;
  delete [] countSendInt;
  delete [] dispSendInt;
  delete [] countRecvInt;
  delete [] dispRecvInt;
}

void freeMPIBuff_dynamic(){
  MPI_Free_mem(sendToRight);
  MPI_Free_mem(sendToLeft);
  MPI_Free_mem(sendToBack);
  MPI_Free_mem(sendToFront);
  MPI_Free_mem(sendToDown);
  MPI_Free_mem(sendToUp);
  MPI_Free_mem(grainNumSendToRight);
  MPI_Free_mem(grainNumSendToLeft);
  MPI_Free_mem(grainNumSendToBack);
  MPI_Free_mem(grainNumSendToFront);
  MPI_Free_mem(grainNumSendToDown);
  MPI_Free_mem(grainNumSendToUp);
  MPI_Free_mem(recvFromLeft);
  MPI_Free_mem(recvFromRight);
  MPI_Free_mem(recvFromFront);
  MPI_Free_mem(recvFromBack);
  MPI_Free_mem(recvFromUp);
  MPI_Free_mem(recvFromDown);
  MPI_Free_mem(grainNumRecvFromRight);
  MPI_Free_mem(grainNumRecvFromLeft);
  MPI_Free_mem(grainNumRecvFromBack);
  MPI_Free_mem(grainNumRecvFromFront);
  MPI_Free_mem(grainNumRecvFromDown);
  MPI_Free_mem(grainNumRecvFromUp);
}

void updateAllPoints(){
  for(int i = 0; i < num_grains; ++i){
    if(belongToThisRank[i] != 0){
      grainsWorld[i].updatePoints();
    }
  }
}
void clearForces(int step){
  for(int i = 0; i < num_grains; ++i){
    if(belongToThisRank[i] != 0){
      grainsWorld[i].changeForce(Vector3d(0.,0.,0.));
      grainsWorld[i].changeMoment(Vector3d(0.,0.,0.));
      grainsWorld[i].changeDOmega(Vector3d(0.,0.,0.));
      if(step % 500 == 0){
         grainsWorld[i].changeVelocity(Vector3d(0.,0.,0.));
         grainsWorld[i].changeOmega(Vector3d(0.,0.,0.));
      }
    }
  }
}
void sort(int & dimX, int & dimY, int & dimZ){
  //so that dimZ >= dimX >= dimY
  int temp;
  if(dimZ < dimY){
    temp = dimZ; dimZ = dimY; dimY = temp;
  }
  if(dimZ < dimX){
    temp = dimZ; dimZ = dimX; dimX = temp;
  }
  if(dimX < dimY){
    temp = dimX; dimX = dimY; dimY = temp;
  }
}
void factorization(int N, int & dimX, int & dimY){
  dimX = 1;
  for(int i = 1; i <= ceil(sqrt(N)); ++i){
    if(N % i == 0) dimX = i;
  }
  dimY = N/dimX;
}
void factorization(int N, int & dimX, int & dimY, int & dimZ){
  dimZ = 1; //it is weird, pow(N, 1./3.) does not work, pow(64,1.0/3)=3.9999
  for(int i = 1; i <= ceil(pow(N, 1.0/3)); ++i){
    if(N % i == 0) dimZ = i;
  }
  factorization(N/dimZ, dimX, dimY);
}
void determine_proc(int num_grains, int num_procs, int & USED, int & dimX, int & dimY, int & dimZ){
  USED = ceil((double)num_grains/1);
  if(USED > num_procs) USED = num_procs;
  int _dimX, _dimY, _dimZ, _USED = USED;

  factorization(USED, dimX, dimY, dimZ);
  while(_USED > 0){
    factorization(_USED, _dimX, _dimY, _dimZ);
    if((double)_USED/(_dimX+_dimY+_dimZ) > (double)USED/(dimX+dimY+dimZ)){
      dimX = _dimX;
      dimY = _dimY;
      dimZ = _dimZ;
    USED = _USED;
    }
    _USED--;
  }
     //USED = 27;
     //dimZ = 3;
     //dimY = 3;
     //dimX = 3;
}
#endif /*MPI_HELPER_FUNCTION_H_*/
