/*
 * MPI_GATHER_FOR_SAVE_H_
 * Created on: June 13, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_GATHER_FOR_SAVE_H_
#define MPI_GATHER_FOR_SAVE_H_
#include "mpi_helper_function.h"

void gather_for_save(int rank, FILE* posfile, FILE* rotfile, string outdir, int step){
  if(rank < USED){
    if(rank > 0){
      MPI_Request request;
      vector<output_data> sendInfo;
      for(int i = 0; i < num_grains; ++i){
        if(belongToThisRank[i] == 1){
          sendInfo.emplace_back(i,grainsWorld[i].getQuat(),grainsWorld[i].getPosition());
        }
      }
      int numOfSendInfo = sendInfo.size();
      MPI_Isend(&numOfSendInfo,1,MPI_INT,0,rank,CART_COMM,&request);
      if(numOfSendInfo != 0){
        MPI_Send(sendInfo.data(),numOfSendInfo,OUTPUT_DATA,0,rank,CART_COMM);
      }
      MPI_Wait(&request, MPI_STATUS_IGNORE);
    }
    if(rank == 0) cout<<"output at step: "<<step<<endl;
    if(rank == 0){
      for(int _rank = 1; _rank < USED; ++_rank){
        MPI_Status r_status;
        int numOfRecvInfo = 0;
        MPI_Recv(&numOfRecvInfo,1,MPI_INT,_rank,_rank,CART_COMM,&r_status);
        if(numOfRecvInfo != 0){
          MPI_Status delivery_status;
          vector<output_data> recvInfo;
          recvInfo.resize(numOfRecvInfo);
          MPI_Recv(recvInfo.data(),numOfRecvInfo,OUTPUT_DATA,_rank,_rank,CART_COMM,&delivery_status);
          for(int i = 0; i < numOfRecvInfo; ++i){
            grainsWorld[recvInfo[i]._id].changeQuat(recvInfo[i]._quat);
            grainsWorld[recvInfo[i]._id].changePosition(recvInfo[i]._position);
          }
        }
      }
      //output
      //cout<<"reach here..."<<endl;
      for (int i = 0; i < num_grains; ++i) {
        const Vector3d position = grainsWorld[i].getPosition();
        const Vector4d rotation = grainsWorld[i].getQuat();
        fprintf(posfile, "%.4f %.4f %.4f\n", position(0), position(1), position(2));
        fprintf(rotfile, "%.4f %.4f %.4f %.4f\n", rotation(0), rotation(1), rotation(2), rotation(3));
      }
      fflush(posfile);
      fflush(rotfile);
    }
  }
}
#endif /*MPI_GATHER_FOR_SAVE_H_*/
