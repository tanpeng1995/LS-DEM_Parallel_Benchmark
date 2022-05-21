#ifndef MPI_MOVE_GRAIN_HPP_
#define MPI_MOVE_GRAIN_HPP_
#include "mpi_helper_function.hpp"

void grainsExchange(int rank, string indir){
  //step 1:
  grain_send_count[LEFT]  = grain_to_be_send[LEFT].size();
  grain_send_count[RIGHT] = grain_to_be_send[RIGHT].size();

  //send to left and receive from right
  MPI_Status status;
  MPI_Sendrecv(&grain_send_count[LEFT],1,MPI_INT,neighbor[LEFT],0,
    &grain_recv_count[RIGHT],1,MPI_INT,neighbor[RIGHT],0,CART_COMM,&status);
  grain_to_be_recv[RIGHT].resize(grain_recv_count[RIGHT]);
  MPI_Sendrecv(grain_to_be_send[LEFT].data(),grain_send_count[LEFT],DATA,neighbor[LEFT],1,
    grain_to_be_recv[RIGHT].data(),grain_recv_count[RIGHT],DATA,neighbor[RIGHT],1,CART_COMM,&status);
  //send to right and receive from left
  MPI_Sendrecv(&grain_send_count[RIGHT],1,MPI_INT,neighbor[RIGHT],2,
    &grain_recv_count[LEFT],1,MPI_INT,neighbor[LEFT],2,CART_COMM,&status);
  grain_to_be_recv[LEFT].resize(grain_recv_count[LEFT]);
  MPI_Sendrecv(grain_to_be_send[RIGHT].data(),grain_send_count[RIGHT],DATA,neighbor[RIGHT],3,
    grain_to_be_recv[LEFT].data(),grain_recv_count[LEFT],DATA,neighbor[LEFT],3,CART_COMM,&status);

  //receive from right
  for(int i = 0; i < grain_recv_count[RIGHT]; ++i){
    int gid = grain_to_be_recv[RIGHT][i]._id;
    int _idx, _idy, _idz, _proc_id;
    find_pos(grain_to_be_recv[RIGHT][i]._position, _idx, _idy, _idz);
    _proc_id = to_proc_id(_idx,_idy,_idz);
    if(_proc_id == rank){
      int _bin_id = to_bin_id(_idx,_idy,_idz);
      bins[_bin_id].insert(gid);

      Vector3d position = grain_to_be_recv[RIGHT][i]._position;
      Vector3d velocity = grain_to_be_recv[RIGHT][i]._velocity;
      Vector3d omega = grain_to_be_recv[RIGHT][i]._omega;
      Vector4d quat(grain_to_be_recv[RIGHT][i]._quat[0], grain_to_be_recv[RIGHT][i]._quat[1],
                    grain_to_be_recv[RIGHT][i]._quat[2], grain_to_be_recv[RIGHT][i]._quat[3]);
      if(grainsWorld[gid].isRealGrain()){
        grainsWorld[gid].changePosition(grain_to_be_recv[RIGHT][i]._position); // update _position
        grainsWorld[gid].changeRotation(grain_to_be_recv[RIGHT][i]._quat); // update _quat, update _rotMatrix, update _pointList
        grainsWorld[gid].changeVelocity(grain_to_be_recv[RIGHT][i]._velocity); // update _velocity
        grainsWorld[gid].changeOmega(grain_to_be_recv[RIGHT][i]._omega); // update _omega, update _omegaGlobal
      }else{
        stringstream morphfile;
        morphfile << indir << "Morphologies/morph_" << gid%target_grains+1 << ".dat";
        grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), grainDensity, gid);
        grainsWorld[gid].setRealGrain(true);
      }
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
      //END ADD
    }else if(_proc_id%(dimX*dimY)/dimX < rank%(dimX*dimY)/dimX){
      grain_to_be_send[FRONT].push_back(grain_to_be_recv[RIGHT][i]);
    }else if(_proc_id%(dimX*dimY)/dimX > rank%(dimX*dimY)/dimX){
      grain_to_be_send[BACK].push_back(grain_to_be_recv[RIGHT][i]);
    }
  }

  //receive from left
  for(int i = 0; i < grain_recv_count[LEFT]; ++i){
    int gid = grain_to_be_recv[LEFT][i]._id;
    int _idx, _idy, _idz, _proc_id;
    find_pos(grain_to_be_recv[LEFT][i]._position, _idx, _idy, _idz);
    _proc_id = to_proc_id(_idx,_idy,_idz);
    if(_proc_id == rank){
      int _bin_id = to_bin_id(_idx,_idy,_idz);
      bins[_bin_id].insert(gid);

      Vector3d position = grain_to_be_recv[LEFT][i]._position;
      Vector3d velocity = grain_to_be_recv[LEFT][i]._velocity;
      Vector3d omega = grain_to_be_recv[LEFT][i]._omega;
      Vector4d quat(grain_to_be_recv[LEFT][i]._quat[0], grain_to_be_recv[LEFT][i]._quat[1],
                    grain_to_be_recv[LEFT][i]._quat[2], grain_to_be_recv[LEFT][i]._quat[3]);
      if(grainsWorld[gid].isRealGrain()){
        grainsWorld[gid].changePosition(grain_to_be_recv[LEFT][i]._position); // update _position
        grainsWorld[gid].changeRotation(grain_to_be_recv[LEFT][i]._quat); // update _quat, update _rotMatrix, update _pointList
        grainsWorld[gid].changeVelocity(grain_to_be_recv[LEFT][i]._velocity); // update _velocity
        grainsWorld[gid].changeOmega(grain_to_be_recv[LEFT][i]._omega); // update _omega, update _omegaGlobal
      }else{
        stringstream morphfile;
        morphfile << indir << "Morphologies/morph_" << gid%target_grains+1 << ".dat";
        grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), grainDensity, gid);
        grainsWorld[gid].setRealGrain(true);
      }
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
      //END ADD
    }else if(_proc_id%(dimX*dimY)/dimX < rank%(dimX*dimY)/dimX){
      grain_to_be_send[FRONT].push_back(grain_to_be_recv[LEFT][i]);
    }else if(_proc_id%(dimX*dimY)/dimX > rank%(dimX*dimY)/dimX){
      grain_to_be_send[BACK].push_back(grain_to_be_recv[LEFT][i]);
    }
  }

  //step 2:
  grain_send_count[FRONT]  = grain_to_be_send[FRONT].size();
  grain_send_count[BACK] = grain_to_be_send[BACK].size();

  //send to front and receive from back
  MPI_Sendrecv(&grain_send_count[FRONT],1,MPI_INT,neighbor[FRONT],0,
    &grain_recv_count[BACK],1,MPI_INT,neighbor[BACK],0,CART_COMM,&status);
  grain_to_be_recv[BACK].resize(grain_recv_count[BACK]);
  MPI_Sendrecv(grain_to_be_send[FRONT].data(),grain_send_count[FRONT],DATA,neighbor[FRONT],1,
    grain_to_be_recv[BACK].data(),grain_recv_count[BACK],DATA,neighbor[BACK],1,CART_COMM,&status);
  //send to back and receive from front
  MPI_Sendrecv(&grain_send_count[BACK],1,MPI_INT,neighbor[BACK],2,
    &grain_recv_count[FRONT],1,MPI_INT,neighbor[FRONT],2,CART_COMM,&status);
  grain_to_be_recv[FRONT].resize(grain_recv_count[FRONT]);
  MPI_Sendrecv(grain_to_be_send[BACK].data(),grain_send_count[BACK],DATA,neighbor[BACK],3,
    grain_to_be_recv[FRONT].data(),grain_recv_count[FRONT],DATA,neighbor[FRONT],3,CART_COMM,&status);

  //receive from back
  for(int i = 0; i < grain_recv_count[BACK]; ++i){
      int gid = grain_to_be_recv[BACK][i]._id;
      int _idx, _idy, _idz, _proc_id;
      find_pos(grain_to_be_recv[BACK][i]._position, _idx, _idy, _idz);
      _proc_id = to_proc_id(_idx,_idy,_idz);
      if(_proc_id == rank){
        int _bin_id = to_bin_id(_idx,_idy,_idz);
        bins[_bin_id].insert(gid);

        Vector3d position = grain_to_be_recv[BACK][i]._position;
        Vector3d velocity = grain_to_be_recv[BACK][i]._velocity;
        Vector3d omega = grain_to_be_recv[BACK][i]._omega;
        Vector4d quat(grain_to_be_recv[BACK][i]._quat[0], grain_to_be_recv[BACK][i]._quat[1],
                      grain_to_be_recv[BACK][i]._quat[2], grain_to_be_recv[BACK][i]._quat[3]);
        if(grainsWorld[gid].isRealGrain()){
          grainsWorld[gid].changePosition(grain_to_be_recv[BACK][i]._position); // update _position
          grainsWorld[gid].changeRotation(grain_to_be_recv[BACK][i]._quat); // update _quat, update _rotMatrix, update _pointList
          grainsWorld[gid].changeVelocity(grain_to_be_recv[BACK][i]._velocity); // update _velocity
          grainsWorld[gid].changeOmega(grain_to_be_recv[BACK][i]._omega); // update _omega, update _omegaGlobal
        }else{
          stringstream morphfile;
          morphfile << indir << "Morphologies/morph_" << gid%target_grains+1 << ".dat";
          grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), grainDensity, gid);
          grainsWorld[gid].setRealGrain(true);
        }
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
        //END ADD
      }else if(_proc_id/(dimX*dimY) < rank/(dimX*dimY)){
        grain_to_be_send[DOWN].push_back(grain_to_be_recv[BACK][i]);
      }else if(_proc_id/(dimX*dimY) > rank/(dimX*dimY)){
        grain_to_be_send[UP].push_back(grain_to_be_recv[BACK][i]);
      }
    }

  //receive from front
  for(int i = 0; i < grain_recv_count[FRONT]; ++i){
      int gid = grain_to_be_recv[FRONT][i]._id;
      int _idx, _idy, _idz, _proc_id;
      find_pos(grain_to_be_recv[FRONT][i]._position, _idx, _idy, _idz);
      _proc_id = to_proc_id(_idx,_idy,_idz);
      if(_proc_id == rank){
        int _bin_id = to_bin_id(_idx,_idy,_idz);
        bins[_bin_id].insert(gid);
        Vector3d position = grain_to_be_recv[FRONT][i]._position;
        Vector3d velocity = grain_to_be_recv[FRONT][i]._velocity;
        Vector3d omega = grain_to_be_recv[FRONT][i]._omega;
        Vector4d quat(grain_to_be_recv[FRONT][i]._quat[0], grain_to_be_recv[FRONT][i]._quat[1],
                      grain_to_be_recv[FRONT][i]._quat[2], grain_to_be_recv[FRONT][i]._quat[3]);
        if(grainsWorld[gid].isRealGrain()){
          grainsWorld[gid].changePosition(grain_to_be_recv[FRONT][i]._position); // update _position
          grainsWorld[gid].changeRotation(grain_to_be_recv[FRONT][i]._quat); // update _quat, update _rotMatrix, update _pointList
          grainsWorld[gid].changeVelocity(grain_to_be_recv[FRONT][i]._velocity); // update _velocity
          grainsWorld[gid].changeOmega(grain_to_be_recv[FRONT][i]._omega); // update _omega, update _omegaGlobal
        }else{
          stringstream morphfile;
          morphfile << indir << "Morphologies/morph_" << gid%target_grains+1 << ".dat";
          grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), grainDensity, gid);
          grainsWorld[gid].setRealGrain(true);
        }
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
        //END ADD
      }else if(_proc_id/(dimX*dimY) < rank/(dimX*dimY)){
        grain_to_be_send[DOWN].push_back(grain_to_be_recv[FRONT][i]);
      }else if(_proc_id/(dimX*dimY) > rank/(dimX*dimY)){
        grain_to_be_send[UP].push_back(grain_to_be_recv[FRONT][i]);
      }
    }

  //step 3
  grain_send_count[DOWN]  = grain_to_be_send[DOWN].size();
  grain_send_count[UP] = grain_to_be_send[UP].size();

  //send to down and receive from up
  MPI_Sendrecv(&grain_send_count[DOWN],1,MPI_INT,neighbor[DOWN],0,
    &grain_recv_count[UP],1,MPI_INT,neighbor[UP],0,CART_COMM,&status);
  grain_to_be_recv[UP].resize(grain_recv_count[UP]);
  MPI_Sendrecv(grain_to_be_send[DOWN].data(),grain_send_count[DOWN],DATA,neighbor[DOWN],1,
    grain_to_be_recv[UP].data(),grain_recv_count[UP],DATA,neighbor[UP],1,CART_COMM,&status);
  //send to up and receive from down
  MPI_Sendrecv(&grain_send_count[UP],1,MPI_INT,neighbor[UP],2,
    &grain_recv_count[DOWN],1,MPI_INT,neighbor[DOWN],2,CART_COMM,&status);
  grain_to_be_recv[DOWN].resize(grain_recv_count[DOWN]);
  MPI_Sendrecv(grain_to_be_send[UP].data(),grain_send_count[UP],DATA,neighbor[UP],3,
    grain_to_be_recv[DOWN].data(),grain_recv_count[DOWN],DATA,neighbor[DOWN],3,CART_COMM,&status);

  //receive from up
  for(int i = 0; i < grain_recv_count[UP]; ++i){
      int gid = grain_to_be_recv[UP][i]._id;
      int _idx, _idy, _idz, _proc_id;
      find_pos(grain_to_be_recv[UP][i]._position, _idx, _idy, _idz);
      _proc_id = to_proc_id(_idx,_idy,_idz);
      if(_proc_id == rank){
        Vector3d position = grain_to_be_recv[UP][i]._position;
        Vector3d velocity = grain_to_be_recv[UP][i]._velocity;
        Vector3d omega = grain_to_be_recv[UP][i]._omega;
        Vector4d quat(grain_to_be_recv[UP][i]._quat[0], grain_to_be_recv[UP][i]._quat[1],
                      grain_to_be_recv[UP][i]._quat[2], grain_to_be_recv[UP][i]._quat[3]);
        if(grainsWorld[gid].isRealGrain()){
          grainsWorld[gid].changePosition(grain_to_be_recv[UP][i]._position); // update _position
          grainsWorld[gid].changeRotation(grain_to_be_recv[UP][i]._quat); // update _quat, update _rotMatrix, update _pointList
          grainsWorld[gid].changeVelocity(grain_to_be_recv[UP][i]._velocity); // update _velocity
          grainsWorld[gid].changeOmega(grain_to_be_recv[UP][i]._omega); // update _omega, update _omegaGlobal
        }else{
          stringstream morphfile;
          morphfile << indir << "Morphologies/morph_" << gid%target_grains+1 << ".dat";
          grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), grainDensity, gid);
          grainsWorld[gid].setRealGrain(true);
        }
        //ADD
        int _bin_id = to_bin_id(_idx,_idy,_idz);
        bins[_bin_id].insert(gid);
        if(firstGrainIdList[_bin_id] == -1){
          firstGrainIdList[_bin_id] = gid;
        }else{
          int j = firstGrainIdList[_bin_id];
          if(gid == 1433){ cout<<"before j is: "<<j<<" bin id is: "<<_bin_id<<endl; }
          while(grainIdList[j] != -1){
            if(gid == 1433){ cout<<"current j is: "<<j<<endl; }
            j = grainIdList[j];
          }
          grainIdList[j] = gid;
          if(gid == 1433){ cout<<"now j becomes: "<<j<<endl; }
        }
        binIdList[gid] = _bin_id;
        belongToThisRank[gid] = 1;
        //END ADD
      }
    }

  //receive from DOWN
  for(int i = 0; i < grain_recv_count[DOWN]; ++i){
      int gid = grain_to_be_recv[DOWN][i]._id;
      int _idx, _idy, _idz, _proc_id;
      find_pos(grain_to_be_recv[DOWN][i]._position, _idx, _idy, _idz);
      _proc_id = to_proc_id(_idx,_idy,_idz);
      if(_proc_id == rank){
        int _bin_id = to_bin_id(_idx,_idy,_idz);
        bins[_bin_id].insert(gid);
        Vector3d position = grain_to_be_recv[DOWN][i]._position;
        Vector3d velocity = grain_to_be_recv[DOWN][i]._velocity;
        Vector3d omega = grain_to_be_recv[DOWN][i]._omega;
        Vector4d quat(grain_to_be_recv[DOWN][i]._quat[0], grain_to_be_recv[DOWN][i]._quat[1],
                      grain_to_be_recv[DOWN][i]._quat[2], grain_to_be_recv[DOWN][i]._quat[3]);
        if(grainsWorld[gid].isRealGrain()){
          grainsWorld[gid].changePosition(grain_to_be_recv[DOWN][i]._position); // update _position
          grainsWorld[gid].changeRotation(grain_to_be_recv[DOWN][i]._quat); // update _quat, update _rotMatrix, update _pointList
          grainsWorld[gid].changeVelocity(grain_to_be_recv[DOWN][i]._velocity); // update _velocity
          grainsWorld[gid].changeOmega(grain_to_be_recv[DOWN][i]._omega); // update _omega, update _omegaGlobal
        }else{
          stringstream morphfile;
          morphfile << indir << "Morphologies/morph_" << gid%target_grains+1 << ".dat";
          grainsWorld[gid] = regenerateGrainFromFile(position, quat, velocity, omega, morphfile.str(), grainDensity, gid);
          grainsWorld[gid].setRealGrain(true);
        }
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
        //END ADD
      }
    }

  //reset;
  for(int i = 0; i < 6; ++i){
    grain_to_be_send[i].clear();
    grain_to_be_send[i].shrink_to_fit();
    grain_to_be_recv[i].clear();
    grain_to_be_recv[i].shrink_to_fit();
    grain_send_count[i] = 0;
    grain_recv_count[i] = 0;
  }
}
#endif /*MPI_MOVE_GRAIN_HPP_*/
