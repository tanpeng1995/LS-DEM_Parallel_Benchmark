/*
 * MPI_BORDER_COMMUNICATION_H_
 * Created on: June 13, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_BORDER_COMMUNICATION_H_
#define MPI_BORDER_COMMUNICATION_H_
#include "mpi_helper_function.h"
#include "mpi_grains.h"
void rightSendData(int * numGrains, data_t (* sendBuff)[SIZE]){
  int _id, _bin_id;
  for(int i = 0; i < binsInDimZ-2; ++i){
    for(int j = 0; j < binsInDimY-2; ++j){
      _id = i*(binsInDimY-2) + j;
      _bin_id = (i+1)*binsInDimX*binsInDimY + (j+1)*binsInDimX + binsInDimX-2;
      numGrains[_id] = bins[_bin_id].size();
      if(numGrains[_id] > SIZE){
        cout<<"SegFault: SIZE OVERFLOW IN MPI_BORDER_COMMUNICATION_H_"<<endl;
        cout<<"POSSIBLE REASON: CUTOFF WAS CHOSEN TOO LARGE, SO THAT #GRAINS IN EACH BINS SHOULD INCREASE."<<endl;
      }
      int k = 0;
      for(std::set<int>::iterator it = bins[_bin_id].begin(); it != bins[_bin_id].end(); ++it){
        int gid = *it;
        sendBuff[_id][k++] = data_t(gid,grainsWorld[gid].getQuat(),
          grainsWorld[gid].getPosition(),grainsWorld[gid].getVelocity(),grainsWorld[gid].getOmega());
      }
    }
  }
}
void leftSendData(int * numGrains, data_t (* sendBuff)[SIZE]){
  int _id, _bin_id;
  for(int i = 0; i < binsInDimZ-2; ++i){
    for(int j = 0; j < binsInDimY-2; ++j){
      _id = i*(binsInDimY-2) + j;
      _bin_id = (i+1)*binsInDimX*binsInDimY + (j+1)*binsInDimX + 1;
      numGrains[_id] = bins[_bin_id].size();
      if(numGrains[_id] > SIZE){
        cout<<"SegFault: SIZE OVERFLOW IN MPI_BORDER_COMMUNICATION_H_"<<endl;
        cout<<"POSSIBLE REASON: CUTOFF WAS CHOSEN TOO LARGE, SO THAT #GRAINS IN EACH BINS SHOULD INCREASE."<<endl;
      }
      int k = 0;
      for(std::set<int>::iterator it = bins[_bin_id].begin(); it != bins[_bin_id].end(); ++it){
        int gid = *it;
        sendBuff[_id][k++] = data_t(gid,grainsWorld[gid].getQuat(),
          grainsWorld[gid].getPosition(),grainsWorld[gid].getVelocity(),grainsWorld[gid].getOmega());
      }
    }
  }
}
void backSendData(int * numGrains, data_t (* sendBuff)[SIZE]){
  int _id, _bin_id;
  for(int i = 0; i < binsInDimZ-2; ++i){
    for(int j = 0; j < binsInDimX; ++j){
      _id = i*binsInDimX + j;
      _bin_id = (i+1)*binsInDimX*binsInDimY + (binsInDimY-2)*binsInDimX + j;
      numGrains[_id] = bins[_bin_id].size();
      if(numGrains[_id] > SIZE){
        cout<<"SegFault: SIZE OVERFLOW IN MPI_BORDER_COMMUNICATION_H_"<<endl;
        cout<<"POSSIBLE REASON: CUTOFF WAS CHOSEN TOO LARGE, SO THAT #GRAINS IN EACH BINS SHOULD INCREASE."<<endl;
      }
      int k = 0;
      for(std::set<int>::iterator it = bins[_bin_id].begin(); it != bins[_bin_id].end(); ++it){
        int gid = *it;
        sendBuff[_id][k++] = data_t(gid,grainsWorld[gid].getQuat(),
          grainsWorld[gid].getPosition(),grainsWorld[gid].getVelocity(),grainsWorld[gid].getOmega());
      }
    }
  }
}
void frontSendData(int * numGrains, data_t (* sendBuff)[SIZE]){
  int _id, _bin_id;
  for(int i = 0; i < binsInDimZ-2; ++i){
    for(int j = 0; j < binsInDimX; ++j){
      _id = i*binsInDimX + j;
      _bin_id = (i+1)*binsInDimX*binsInDimY + binsInDimX + j;
      numGrains[_id] = bins[_bin_id].size();
      if(numGrains[_id] > SIZE){
        cout<<"SegFault: SIZE OVERFLOW IN MPI_BORDER_COMMUNICATION_H_"<<endl;
        cout<<"POSSIBLE REASON: CUTOFF WAS CHOSEN TOO LARGE, SO THAT #GRAINS IN EACH BINS SHOULD INCREASE."<<endl;
      }
      int k = 0;
      for(std::set<int>::iterator it = bins[_bin_id].begin(); it != bins[_bin_id].end(); ++it){
        int gid = *it;
        sendBuff[_id][k++] = data_t(gid,grainsWorld[gid].getQuat(),
          grainsWorld[gid].getPosition(),grainsWorld[gid].getVelocity(),grainsWorld[gid].getOmega());
      }
    }
  }
}
void downSendData(int * numGrains, data_t (* sendBuff)[SIZE]){
  int _id, _bin_id;
  for(int i = 0; i < binsInDimY; ++i){
    for(int j = 0; j < binsInDimX; ++j){
      _id = i*binsInDimX + j;
      _bin_id = binsInDimX*binsInDimY + i*binsInDimX + j;
      numGrains[_id] = bins[_bin_id].size();
      if(numGrains[_id] > SIZE){
        cout<<"SegFault: SIZE OVERFLOW IN MPI_BORDER_COMMUNICATION_H_"<<endl;
        cout<<"POSSIBLE REASON: CUTOFF WAS CHOSEN TOO LARGE, SO THAT #GRAINS IN EACH BINS SHOULD INCREASE."<<endl;
      }
      int k = 0;
      for(std::set<int>::iterator it = bins[_bin_id].begin(); it != bins[_bin_id].end(); ++it){
        int gid = *it;
        sendBuff[_id][k++] = data_t(gid,grainsWorld[gid].getQuat(),
          grainsWorld[gid].getPosition(),grainsWorld[gid].getVelocity(),grainsWorld[gid].getOmega());
      }
    }
  }
}
void upSendData(int * numGrains, data_t (* sendBuff)[SIZE]){
  int _id, _bin_id;
  for(int i = 0; i < binsInDimY; ++i){
    for(int j = 0; j < binsInDimX; ++j){
      _id = i*binsInDimX + j;
      _bin_id = (binsInDimZ-2)*binsInDimX*binsInDimY + i*binsInDimX + j;
      numGrains[_id] = bins[_bin_id].size();
      if(numGrains[_id] > SIZE){
        cout<<"SegFault: SIZE OVERFLOW IN MPI_BORDER_COMMUNICATION_H_"<<endl;
        cout<<"POSSIBLE REASON: CUTOFF WAS CHOSEN TOO LARGE, SO THAT #GRAINS IN EACH BINS SHOULD INCREASE."<<endl;
      }
      int k = 0;
      for(std::set<int>::iterator it = bins[_bin_id].begin(); it != bins[_bin_id].end(); ++it){
        int gid = *it;
        sendBuff[_id][k++] = data_t(gid,grainsWorld[gid].getQuat(),
          grainsWorld[gid].getPosition(),grainsWorld[gid].getVelocity(),grainsWorld[gid].getOmega());
      }
    }
  }
}
void updateRightLeftSide(string indir){
  MPI_Status status;
  rightSendData(grainNumSendToRight,sendToRight);
  leftSendData(grainNumSendToLeft,sendToLeft);
  //send to right and recv from left
  MPI_Sendrecv(grainNumSendToRight,(binsInDimY-2)*(binsInDimZ-2),MPI_INT,neighbor[RIGHT],0,
    grainNumRecvFromLeft,(binsInDimY-2)*(binsInDimZ-2),MPI_INT,neighbor[LEFT],0,CART_COMM,&status);
  MPI_Sendrecv(sendToRight,(binsInDimY-2)*(binsInDimZ-2),BIN,neighbor[RIGHT],1,
    recvFromLeft,(binsInDimY-2)*(binsInDimZ-2),BIN,neighbor[LEFT],1,CART_COMM,&status);
  //send to left and recv from right
  MPI_Sendrecv(grainNumSendToLeft,(binsInDimY-2)*(binsInDimZ-2),MPI_INT,neighbor[LEFT],2,
    grainNumRecvFromRight,(binsInDimY-2)*(binsInDimZ-2),MPI_INT,neighbor[RIGHT],2,CART_COMM,&status);
  MPI_Sendrecv(sendToLeft,(binsInDimY-2)*(binsInDimZ-2),BIN,neighbor[LEFT],3,
    recvFromRight,(binsInDimY-2)*(binsInDimZ-2),BIN,neighbor[RIGHT],3,CART_COMM,&status);

  int _id, _bin_id, gid;
  //update left side
  for(int i = 0; i < binsInDimZ-2; ++i){
    for(int j = 0; j < binsInDimY-2; ++j){
      _id = i*(binsInDimY-2) + j;
      _bin_id = (i+1)*binsInDimX*binsInDimY + (j+1)*binsInDimX;
      bins[_bin_id].clear();
      firstGrainIdList[_bin_id] = -1;
      for(int k = 0; k < grainNumRecvFromLeft[_id]; ++k){
        gid = recvFromLeft[_id][k]._id;
        grainIdList[gid] = -1;
        belongToThisRank[gid] = 2;
        //grainsWorld[gid].changePosition(recvFromLeft[_id][k]._position);
        //grainsWorld[gid].changeRotation(recvFromLeft[_id][k]._quat);
        //grainsWorld[gid].changeVelocity(recvFromLeft[_id][k]._velocity);
        //grainsWorld[gid].changeOmega(recvFromLeft[_id][k]._omega);
        Vector3d position = recvFromLeft[_id][k]._position;
        Vector4d quat(recvFromLeft[_id][k]._quat[0], recvFromLeft[_id][k]._quat[1],
                      recvFromLeft[_id][k]._quat[2], recvFromLeft[_id][k]._quat[3]);
        Vector3d velocity = recvFromLeft[_id][k]._velocity;
        Vector3d omega = recvFromLeft[_id][k]._omega;
        if(grainsWorld[gid].getRealGrain()){
          grainsWorld[gid].changePosition(recvFromLeft[_id][k]._position);
          grainsWorld[gid].changeRotation(recvFromLeft[_id][k]._quat);
          grainsWorld[gid].changeVelocity(recvFromLeft[_id][k]._velocity);
          grainsWorld[gid].changeOmega(recvFromLeft[_id][k]._omega);
        }else{
          //cout<<"read ghost grain: "<<gid<<endl;
          stringstream morphfile;
          morphfile << indir << "Morphologies/morph_" << gid+1 << ".dat";
          grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), gid);
          grainsWorld[gid].changeDensity(grainDensity);
          grainsWorld[gid].setRealGrain(true);
        }

        if(firstGrainIdList[_bin_id] == -1){ // need to clear bin & part list
          firstGrainIdList[_bin_id] = gid;
        }else{
          int j = firstGrainIdList[_bin_id];
          while(grainIdList[j] != -1) j = grainIdList[j];
          grainIdList[j] = gid;
        }
        binIdList[gid] = _bin_id;
        bins[_bin_id].insert(gid);
      }
    }
  }
  for(int i = 0; i < (binsInDimY-2)*(binsInDimZ-2); ++i){
    grainNumSendToRight[i] = 0;
    grainNumRecvFromLeft[i] = 0;
    for(int j = 0; j < SIZE; ++j){
      sendToRight[i][j] = data_t();
      recvFromLeft[i][j] = data_t();
    }
  }
  //update right side
  for(int i = 0; i < binsInDimZ-2; ++i){
    for(int j = 0; j < binsInDimY-2; ++j){
      _id = i*(binsInDimY-2) + j;
      _bin_id = (i+1)*binsInDimX*binsInDimY + (j+1)*binsInDimX + binsInDimX-1;
      bins[_bin_id].clear();
      firstGrainIdList[_bin_id] = -1;
      for(int k = 0; k < grainNumRecvFromRight[_id]; ++k){
        gid = recvFromRight[_id][k]._id;
        //bins[_bin_id].insert(gid);
        grainIdList[gid] = -1;
        belongToThisRank[gid] = 2;
        //grainsWorld[gid].changePosition(recvFromRight[_id][k]._position);
        //grainsWorld[gid].changeRotation(recvFromRight[_id][k]._quat);
        //grainsWorld[gid].changeVelocity(recvFromRight[_id][k]._velocity);
        //grainsWorld[gid].changeOmega(recvFromRight[_id][k]._omega);

        Vector3d position = recvFromRight[_id][k]._position;
        Vector4d quat(recvFromRight[_id][k]._quat[0], recvFromRight[_id][k]._quat[1],
                      recvFromRight[_id][k]._quat[2], recvFromRight[_id][k]._quat[3]);
        Vector3d velocity = recvFromRight[_id][k]._velocity;
        Vector3d omega = recvFromRight[_id][k]._omega;
        if(grainsWorld[gid].getRealGrain()){
          grainsWorld[gid].changePosition(recvFromRight[_id][k]._position);
          grainsWorld[gid].changeRotation(recvFromRight[_id][k]._quat);
          grainsWorld[gid].changeVelocity(recvFromRight[_id][k]._velocity);
          grainsWorld[gid].changeOmega(recvFromRight[_id][k]._omega);
        }else{
          //cout<<"read ghost grain: "<<gid<<endl;
          stringstream morphfile;
          morphfile << indir << "Morphologies/morph_" << gid+1 << ".dat";
          grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), gid);
          grainsWorld[gid].changeDensity(grainDensity);
          grainsWorld[gid].setRealGrain(true);
        }

        if(firstGrainIdList[_bin_id] == -1){ // need to clear bin & part list
          firstGrainIdList[_bin_id] = gid;
        }else{
          int j = firstGrainIdList[_bin_id];
          while(grainIdList[j] != -1) j = grainIdList[j];
          grainIdList[j] = gid;
        }
        binIdList[gid] = _bin_id;
        bins[_bin_id].insert(gid);
      }
    }
  }
  for(int i = 0; i < (binsInDimY-2)*(binsInDimZ-2); ++i){
    grainNumSendToLeft[i] = 0;
    grainNumRecvFromRight[i] = 0;
    for(int j = 0; j < SIZE; ++j){
      sendToLeft[i][j] = data_t();
      recvFromRight[i][j] = data_t();
    }
  }
}
void updateBackFrontSide(string indir){
  MPI_Status status;
  backSendData(grainNumSendToBack,sendToBack);
  frontSendData(grainNumSendToFront,sendToFront);
  //send to back and recv from front
  MPI_Sendrecv(grainNumSendToBack,binsInDimX*(binsInDimZ-2),MPI_INT,neighbor[BACK],0,
    grainNumRecvFromFront,binsInDimX*(binsInDimZ-2),MPI_INT,neighbor[FRONT],0,CART_COMM,&status);
  MPI_Sendrecv(sendToBack,binsInDimX*(binsInDimZ-2),BIN,neighbor[BACK],1,
    recvFromFront,binsInDimX*(binsInDimZ-2),BIN,neighbor[FRONT],1,CART_COMM,&status);
  //send to front and recv from back
  MPI_Sendrecv(grainNumSendToFront,binsInDimX*(binsInDimZ-2),MPI_INT,neighbor[FRONT],2,
    grainNumRecvFromBack,binsInDimX*(binsInDimZ-2),MPI_INT,neighbor[BACK],2,CART_COMM,&status);
  MPI_Sendrecv(sendToFront,binsInDimX*(binsInDimZ-2),BIN,neighbor[FRONT],3,
    recvFromBack,binsInDimX*(binsInDimZ-2),BIN,neighbor[BACK],3,CART_COMM,&status);

  int _id, _bin_id, gid;
  //update front side
  for(int i = 0; i < binsInDimZ-2; ++i){
    for(int j = 0; j < binsInDimX; ++j){
      _id = i*binsInDimX + j;
      _bin_id = (i+1)*binsInDimX*binsInDimY + j;
      bins[_bin_id].clear();
      firstGrainIdList[_bin_id] = -1;
      for(int k = 0; k < grainNumRecvFromFront[_id]; ++k){
        gid = recvFromFront[_id][k]._id;
        //bins[_bin_id].insert(gid);
        grainIdList[gid] = -1;
        belongToThisRank[gid] = 2;
        //grainsWorld[gid].changePosition(recvFromFront[_id][k]._position);
        //grainsWorld[gid].changeRotation(recvFromFront[_id][k]._quat);
        //grainsWorld[gid].changeVelocity(recvFromFront[_id][k]._velocity);
        //grainsWorld[gid].changeOmega(recvFromFront[_id][k]._omega);

        Vector3d position = recvFromFront[_id][k]._position;
        Vector4d quat(recvFromFront[_id][k]._quat[0], recvFromFront[_id][k]._quat[1],
                      recvFromFront[_id][k]._quat[2], recvFromFront[_id][k]._quat[3]);
        Vector3d velocity = recvFromFront[_id][k]._velocity;
        Vector3d omega = recvFromFront[_id][k]._omega;
        if(grainsWorld[gid].getRealGrain()){
          grainsWorld[gid].changePosition(recvFromFront[_id][k]._position);
          grainsWorld[gid].changeRotation(recvFromFront[_id][k]._quat);
          grainsWorld[gid].changeVelocity(recvFromFront[_id][k]._velocity);
          grainsWorld[gid].changeOmega(recvFromFront[_id][k]._omega);
        }else{
          //cout<<"read ghost grain: "<<gid<<endl;
          stringstream morphfile;
          morphfile << indir << "Morphologies/morph_" << gid+1 << ".dat";
          grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), gid);
          grainsWorld[gid].changeDensity(grainDensity);
          grainsWorld[gid].setRealGrain(true);
        }

        if(firstGrainIdList[_bin_id] == -1){ // need to clear bin & part list
          firstGrainIdList[_bin_id] = gid;
        }else{
          int j = firstGrainIdList[_bin_id];
          while(grainIdList[j] != -1) j = grainIdList[j];
          grainIdList[j] = gid;
        }
        binIdList[gid] = _bin_id;
        bins[_bin_id].insert(gid);
      }
    }
  }
  for(int i = 0; i < binsInDimX*(binsInDimZ-2); ++i){
    grainNumSendToBack[i] = 0;
    grainNumRecvFromFront[i] = 0;
    for(int j = 0; j < SIZE; ++j){
      sendToBack[i][j] = data_t();
      recvFromFront[i][j] = data_t();
    }
  }
  //update back side
  for(int i = 0; i < binsInDimZ-2; ++i){
    for(int j = 0; j < binsInDimX; ++j){
      _id = i*binsInDimX + j;
      _bin_id = (i+1)*binsInDimX*binsInDimY + (binsInDimY-1)*binsInDimX + j;
      bins[_bin_id].clear();
      firstGrainIdList[_bin_id] = -1;
      for(int k = 0; k < grainNumRecvFromBack[_id]; ++k){
        gid = recvFromBack[_id][k]._id;
        //bins[_bin_id].insert(gid);
        grainIdList[gid] = -1;
        belongToThisRank[gid] = 2;
        //grainsWorld[gid].changePosition(recvFromBack[_id][k]._position);
        //grainsWorld[gid].changeRotation(recvFromBack[_id][k]._quat);
        //grainsWorld[gid].changeVelocity(recvFromBack[_id][k]._velocity);
        //grainsWorld[gid].changeOmega(recvFromBack[_id][k]._omega);

        Vector3d position = recvFromBack[_id][k]._position;
        Vector4d quat(recvFromBack[_id][k]._quat[0], recvFromBack[_id][k]._quat[1],
                      recvFromBack[_id][k]._quat[2], recvFromBack[_id][k]._quat[3]);
        Vector3d velocity = recvFromBack[_id][k]._velocity;
        Vector3d omega = recvFromBack[_id][k]._omega;
        if(grainsWorld[gid].getRealGrain()){
          grainsWorld[gid].changePosition(recvFromBack[_id][k]._position);
          grainsWorld[gid].changeRotation(recvFromBack[_id][k]._quat);
          grainsWorld[gid].changeVelocity(recvFromBack[_id][k]._velocity);
          grainsWorld[gid].changeOmega(recvFromBack[_id][k]._omega);
        }else{
          //cout<<"read ghost grain: "<<gid<<endl;
          stringstream morphfile;
          morphfile << indir << "Morphologies/morph_" << gid+1 << ".dat";
          grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), gid);
          grainsWorld[gid].changeDensity(grainDensity);
          grainsWorld[gid].setRealGrain(true);
        }

        if(firstGrainIdList[_bin_id] == -1){ // need to clear bin & part list
          firstGrainIdList[_bin_id] = gid;
        }else{
          int j = firstGrainIdList[_bin_id];
          while(grainIdList[j] != -1) j = grainIdList[j];
          grainIdList[j] = gid;
        }
        binIdList[gid] = _bin_id;
        bins[_bin_id].insert(gid);
      }
    }
  }
  for(int i = 0; i < binsInDimX*(binsInDimZ-2); ++i){
    grainNumSendToFront[i] = 0;
    grainNumRecvFromBack[i] = 0;
    for(int j = 0; j < SIZE; ++j){
      sendToFront[i][j] = data_t();
      recvFromBack[i][j] = data_t();
    }
  }
}
void updateDownUpSide(string indir){
  MPI_Status status;
  upSendData(grainNumSendToUp,sendToUp);
  downSendData(grainNumSendToDown,sendToDown);
  //send to up and recv from down
  MPI_Sendrecv(grainNumSendToUp,binsInDimX*binsInDimY,MPI_INT,neighbor[UP],0,
    grainNumRecvFromDown,binsInDimX*binsInDimY,MPI_INT,neighbor[DOWN],0,CART_COMM,&status);
  MPI_Sendrecv(sendToUp,binsInDimX*binsInDimY,BIN,neighbor[UP],1,
    recvFromDown,binsInDimX*binsInDimY,BIN,neighbor[DOWN],1,CART_COMM,&status);
  //send to down and recv from up
  MPI_Sendrecv(grainNumSendToDown,binsInDimX*binsInDimY,MPI_INT,neighbor[DOWN],2,
    grainNumRecvFromUp,binsInDimX*binsInDimY,MPI_INT,neighbor[UP],2,CART_COMM,&status);
  MPI_Sendrecv(sendToDown,binsInDimX*binsInDimY,BIN,neighbor[DOWN],3,
    recvFromUp,binsInDimX*binsInDimY,BIN,neighbor[UP],3,CART_COMM,&status);

  int _id, _bin_id, gid;
  //update down side
  for(int i = 0; i < binsInDimY; ++i){
    for(int j = 0; j < binsInDimX; ++j){
      _id = i*binsInDimX + j;
      _bin_id = i*binsInDimX + j;
      bins[_bin_id].clear();
      firstGrainIdList[_bin_id] = -1;
      for(int k = 0; k < grainNumRecvFromDown[_id]; ++k){
        gid = recvFromDown[_id][k]._id;
        //bins[_bin_id].insert(gid);
        grainIdList[gid] = -1;
        belongToThisRank[gid] = 2;
        //grainsWorld[gid].changePosition(recvFromDown[_id][k]._position);
        //grainsWorld[gid].changeRotation(recvFromDown[_id][k]._quat);
        //grainsWorld[gid].changeVelocity(recvFromDown[_id][k]._velocity);
        //grainsWorld[gid].changeOmega(recvFromDown[_id][k]._omega);

        Vector3d position = recvFromDown[_id][k]._position;
        Vector4d quat(recvFromDown[_id][k]._quat[0], recvFromDown[_id][k]._quat[1],
                      recvFromDown[_id][k]._quat[2], recvFromDown[_id][k]._quat[3]);
        Vector3d velocity = recvFromDown[_id][k]._velocity;
        Vector3d omega = recvFromDown[_id][k]._omega;
        if(grainsWorld[gid].getRealGrain()){
          grainsWorld[gid].changePosition(recvFromDown[_id][k]._position);
          grainsWorld[gid].changeRotation(recvFromDown[_id][k]._quat);
          grainsWorld[gid].changeVelocity(recvFromDown[_id][k]._velocity);
          grainsWorld[gid].changeOmega(recvFromDown[_id][k]._omega);
        }else{
          //cout<<"read ghost grain: "<<gid<<endl;
          stringstream morphfile;
          morphfile << indir << "Morphologies/morph_" << gid+1 << ".dat";
          grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), gid);
          grainsWorld[gid].changeDensity(grainDensity);
          grainsWorld[gid].setRealGrain(true);
        }

        if(firstGrainIdList[_bin_id] == -1){ // need to clear bin & part list
          firstGrainIdList[_bin_id] = gid;
        }else{
          int j = firstGrainIdList[_bin_id];
          while(grainIdList[j] != -1) j = grainIdList[j];
          grainIdList[j] = gid;
        }
        binIdList[gid] = _bin_id;
        bins[_bin_id].insert(gid);
      }
    }
  }
  for(int i = 0; i < binsInDimX*binsInDimY; ++i){
    grainNumSendToUp[i] = 0;
    grainNumRecvFromDown[i] = 0;
    for(int j = 0; j < SIZE; ++j){
      sendToUp[i][j] = data_t();
      recvFromDown[i][j] = data_t();
    }
  }
  //update up side
  for(int i = 0; i < binsInDimY; ++i){
    for(int j = 0; j < binsInDimX; ++j){
      _id = i*binsInDimX + j;
      _bin_id = (binsInDimZ-1)*binsInDimX*binsInDimY + i*binsInDimX + j;
      bins[_bin_id].clear();
      firstGrainIdList[_bin_id] = -1;
      for(int k = 0; k < grainNumRecvFromUp[_id]; ++k){
        gid = recvFromUp[_id][k]._id;
        //bins[_bin_id].insert(gid);
        grainIdList[gid] = -1;
        belongToThisRank[gid] = 2;
        //grainsWorld[gid].changePosition(recvFromUp[_id][k]._position);
        //grainsWorld[gid].changeRotation(recvFromUp[_id][k]._quat);
        //grainsWorld[gid].changeVelocity(recvFromUp[_id][k]._velocity);
        //grainsWorld[gid].changeOmega(recvFromUp[_id][k]._omega);

        Vector3d position = recvFromUp[_id][k]._position;
        Vector4d quat(recvFromUp[_id][k]._quat[0], recvFromUp[_id][k]._quat[1],
                      recvFromUp[_id][k]._quat[2], recvFromUp[_id][k]._quat[3]);
        Vector3d velocity = recvFromUp[_id][k]._velocity;
        Vector3d omega = recvFromUp[_id][k]._omega;
        if(grainsWorld[gid].getRealGrain()){
          grainsWorld[gid].changePosition(recvFromUp[_id][k]._position);
          grainsWorld[gid].changeRotation(recvFromUp[_id][k]._quat);
          grainsWorld[gid].changeVelocity(recvFromUp[_id][k]._velocity);
          grainsWorld[gid].changeOmega(recvFromUp[_id][k]._omega);
        }else{
          //cout<<"read ghost grain: "<<gid<<endl;
          stringstream morphfile;
          morphfile << indir << "Morphologies/morph_" << gid+1 << ".dat";
          grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), gid);
          grainsWorld[gid].changeDensity(grainDensity);
          grainsWorld[gid].setRealGrain(true);
        }

        if(firstGrainIdList[_bin_id] == -1){ // need to clear bin & part list
          firstGrainIdList[_bin_id] = gid;
        }else{
          int j = firstGrainIdList[_bin_id];
          while(grainIdList[j] != -1) j = grainIdList[j];
          grainIdList[j] = gid;
        }
        binIdList[gid] = _bin_id;
        bins[_bin_id].insert(gid);
      }
    }
  }
  for(int i = 0; i < binsInDimX*binsInDimY; ++i){
    grainNumSendToDown[i] = 0;
    grainNumRecvFromUp[i] = 0;
    for(int j = 0; j < SIZE; ++j){
      sendToDown[i][j] = data_t();
      recvFromUp[i][j] = data_t();
    }
  }
}
#endif /*MPI_BORDER_COMMUNICATION_H_*/
