/*
 * MPI_SOLVER_H_
 * Created on: June 21, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_SOLVER_H_
#define MPI_SOLVER_H_
#include "mpi_border_communication.h"
#include "mpi_particle_migration.h"
#include "mpi_ballWalls_contact.h"
#include "tuple"

const int T = 5;

void simulate_one_step(int rank, double dt, int step, string indir, string outdir, FILE * volumeCNFile, double & T_force, double & T_border, double & T_migrate, double & T_comm){
  if(rank <= USED){
    /* change loading stage */
    if(std::find(loadingStages.begin(), loadingStages.end(), step) != loadingStages.end()){
      dHeight *= -1;
      if(rank == 0){
        cout << "################## CHANGE LOADING STAGE ###################" << endl;
        cout << "changed loadingStages at step: "<<step<<" current dHeight is: "<<dHeight<<endl;
        cout << "################## CHANGE LOADING STAGE ###################" << endl<<endl;
      }
    }
    /*########################################################################*/
    /*############ BORDER COMMUNICATION AND PREPARE BOUNDARY DATA ############*/
    /*########################################################################*/
    stressVoigt << 0.,0.,0.,0.,0.,0.;
    auto start_time1 = std::chrono::steady_clock::now();
    if(rank < USED){
      if(USED > 1){ updateRightLeftSide(indir); updateBackFrontSide(indir); updateDownUpSide(indir); }
      if(step % T == 0){ collectBoundaryGrains(boundaryDataSend, CUTOFF, initRadius); }
      clearForces(); updateAllPoints();
    }
    auto end_time1 = std::chrono::steady_clock::now(); std::chrono::duration<double> diff1 = end_time1 - start_time1; double seconds1 = diff1.count();
    T_border += seconds1;
    /*########################################################################*/
    /*############# EXCHANGE FOR MEMBRANE INTERACTION (SEND->) ###############*/
    /*########################################################################*/
    auto start_time4 = std::chrono::steady_clock::now();
    int send, recv;
    if(step % T == 0){
      if(rank < USED){
        MPI_Request request; send = boundaryDataSend.size();
        MPI_Isend(&send, 1 , MPI_INT, USED, rank, BOUNDARY_COMM, &request);
        if(send != 0){ MPI_Send(boundaryDataSend.data(), send, BOUNDARY_SEND, USED, rank, BOUNDARY_COMM); }
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }else{
        for(int i = 0; i < USED; ++i){
          MPI_Status s_status; MPI_Recv(&recv, 1, MPI_INT, i, i, BOUNDARY_COMM, &s_status);
          if(recv != 0){
            MPI_Status delivery_status; totalBoundaryDataSend[i].resize(recv);
            MPI_Recv(totalBoundaryDataSend[i].data(), recv, BOUNDARY_SEND, i, i, BOUNDARY_COMM, &delivery_status);
          }
        }
      }
      /* end */
    }
    auto end_time4 = std::chrono::steady_clock::now(); std::chrono::duration<double> diff4 = end_time4 - start_time4; double seconds4 = diff4.count();
    T_comm += seconds4;

    /*########################################################################*/
    /*########################## FORCE RESOLUTION ############################*/
    /*########################################################################*/
    vector<Vector6d> posChain;
    vector<Vector3d> forceChain;
    auto start_time2 = std::chrono::steady_clock::now();
    int j, _idx, _idy, _idz, _bin_id, _proc_id, _old_bin_id;
    if(rank < USED){
      #pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads) private(j, _idx, _idy, _idz, _bin_id, _proc_id, _old_bin_id) shared(grainsWorld, belongToThisRank, grainIdList, wallCap, wallPlane, firstGrainIdList, forceChain, posChain)
      for(int i = 0; i < num_grains; ++i){
        /*
           --|----|----|----|--
             |  6 |  7 |  8 |
           --|----|----|----|--
             |  3 |  4 |  5 |
           --|----|----|----|--
             |  0 |  1 |  2 |
           --|----|----|----|--
        */
        int binNum = binIdList[i];
        int binArray[13] = {binNum-1, binNum+binsInDimX-1, binNum+binsInDimX, binNum+binsInDimX+1,
          binNum+binsInDimX*binsInDimY-binsInDimX-1, binNum+binsInDimX*binsInDimY-binsInDimX,
          binNum+binsInDimX*binsInDimY-binsInDimX+1, binNum+binsInDimX*binsInDimY-1,
          binNum+binsInDimX*binsInDimY, binNum+binsInDimX*binsInDimY+1,
          binNum+binsInDimX*binsInDimY+binsInDimX-1, binNum+binsInDimX*binsInDimY+binsInDimX,
          binNum+binsInDimX*binsInDimY+binsInDimX+1};

        if(belongToThisRank[i]==1){
          j = grainIdList[i];
          //inside the same bins
          while(j != -1){
            if(grainsWorld[i]->bCircleCheck(grainsWorld[j])){
              if(i < j){ grainsWorld[i]->findInterparticleForceMoment(grainsWorld[j],dt,stressVoigt,posChain,forceChain,true); }
              else{ grainsWorld[j]->findInterparticleForceMoment(grainsWorld[i],dt,stressVoigt,posChain,forceChain,true); }
            }
            j = grainIdList[j];
          }
          //consider 3, 6, 7, 8 in the same level as the bins of interest and 0~9 at upper level
          for (int k = 0; k < 13; ++k) {
            j = firstGrainIdList[binArray[k]];
            while(j != -1){
              if(grainsWorld[i]->bCircleCheck(grainsWorld[j])){
                if(i < j){ grainsWorld[i]->findInterparticleForceMoment(grainsWorld[j],dt,stressVoigt,posChain,forceChain,true); }
                else{ grainsWorld[j]->findInterparticleForceMoment(grainsWorld[i],dt,stressVoigt,posChain,forceChain,true); }
              }
              j = grainIdList[j];
            }
          }
          /*####################### CONTACT WITH CAP #########################*/
          if(wallCap->bCircleCheck(grainsWorld[i])){
            wallCap->findWallForceMoment(grainsWorld[i],dt,stressVoigt);
          }
          /*###################### CONTACT WITH BOTTOM #######################*/
          if(wallPlane->bCircleCheck(grainsWorld[i])){
            wallPlane->findWallForceMoment(grainsWorld[i],dt,stressVoigt);
          }

        } else if(belongToThisRank[i] == 2) {
          for (int k = 0; k < 13; ++k) {
            find_pos(binArray[k], _idx, _idy, _idz);
            if(_idz < binsInDimZ-1 && _idz > 0 && _idy < binsInDimY-1 && _idy > 0 && _idx < binsInDimX-1 && _idx > 0) {
              j = firstGrainIdList[binArray[k]];
              while(j != -1){
                if(belongToThisRank[j] == 1 && grainsWorld[i]->bCircleCheck(grainsWorld[j])){
                  if(i < j){ grainsWorld[i]->findInterparticleForceMoment(grainsWorld[j],dt,stressVoigt,posChain,forceChain,false); }
                  else{ grainsWorld[j]->findInterparticleForceMoment(grainsWorld[i],dt,stressVoigt,posChain,forceChain,false); }
                }
                j = grainIdList[j];
              }
            }
          }
        }
      }
    }else{
      /*###################### COMPUTE MEMBRANE INTERACTION IN MEMBRANE PROCESSOR #######################*/
      /* COMPUTE MEMBRANE INTERACTION WITH GRAINS */

      #pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads) shared(grainsWorld, totalBoundaryDataSend, wallBalls, externalWall, totalBoundaryDataRecv)
      for(int i = 0; i < USED; ++i){
        if(step % T == 0){
          for(int j = 0; j < totalBoundaryDataSend[i].size(); ++j){
            int gid = totalBoundaryDataSend[i][j]._id;
            Vector3d force(0.,0.,0.);
            Vector3d moment(0.,0.,0.);
            Vector3d position = totalBoundaryDataSend[i][j]._position;
            Vector4d quat(totalBoundaryDataSend[i][j]._quat[0], totalBoundaryDataSend[i][j]._quat[1],
                          totalBoundaryDataSend[i][j]._quat[2], totalBoundaryDataSend[i][j]._quat[3]);
            if(!grainsWorld[gid]->getRealGrain()){
              stringstream morphfile;
              //cout<<"read grain "<<gid<<" from database to membrane processor."<<endl;
              morphfile << indir << "Morphologies/morph_" << gid+1 << ".dat";
              grainsWorld[gid] = regenerateGrainFromFile(position, quat, Vector3d(0.,0.,0.), Vector3d(0.,0.,0.), morphfile.str(), gid);
              grainsWorld[gid]->changeDensity(grainDensity);
              grainsWorld[gid]->setRealGrain(true);
              if(massScaling){
                if(grainsWorld[gid]->getMass() < minMass){
                  grainsWorld[gid]->changeMass(minMass);
                }
              }
            }else{
              grainsWorld[gid]->changePosition(totalBoundaryDataSend[i][j]._position);
              grainsWorld[gid]->changeRotation(totalBoundaryDataSend[i][j]._quat);
            }
            grainsWorld[gid]->updatePoints();
            //if(wallBalls->bCircleCheck(grainsWorld[gid])){
              wallBalls->findWallForceMoment(grainsWorld[gid],force,moment,stressVoigt);
            //}
            //if(externalWall->bCircleCheck(grainsWorld[gid])){
              externalWall->findWallForceMoment(grainsWorld[gid],force,moment,stressVoigt);
            //}
            totalBoundaryDataRecv[i].emplace_back(gid, force, moment);
          }
        }
        totalBoundaryDataSend[i].clear(); totalBoundaryDataSend.shrink_to_fit();
      }
      /* COMPUTE MEMBRANE INTERACTION WITH CAP */
      wallBalls->computeInternalForceOnTopLayer(wallCap, dt);
      externalWall->computeInternalForceOnTopLayer(wallCap, dt);

    }
    auto end_time2 = std::chrono::steady_clock::now(); std::chrono::duration<double> diff2 = end_time2 - start_time2; double seconds2 = diff2.count();
    T_force += seconds2;

    /*########################################################################*/
    /*############### FIND CAP AND PLANE FORCES WITH GRAINS ##################*/
    /*########################################################################*/
    MPI_Barrier(BOUNDARY_COMM);
    if(rank < USED){
      Vector3d wallCapForce   = wallCap->getForceNonConst();
      Vector3d wallPlaneForce = wallPlane->getForceNonConst();
      MPI_Allreduce(MPI_IN_PLACE, wallCapForce.data(), 3, MPI_DOUBLE, MPI_SUM, CART_COMM);
      MPI_Allreduce(MPI_IN_PLACE, wallPlaneForce.data(), 3, MPI_DOUBLE, MPI_SUM, CART_COMM);
      nowCapPressure = abs( wallCapForce.dot(wallCap->getNormalGlobal()) )/ wallCap->getArea();
      planePressure  = abs( wallPlaneForce.dot(wallPlane->getNormalGlobal()) )/ wallCap->getArea();
      if(rank == 0){ cout<<"current cap pressure is: "<<nowCapPressure/scalingFactor<<" last cap pressure is: "<<prevCapPressure/scalingFactor<<" planePressure is: "<<planePressure/scalingFactor<<endl; }
      prevCapPressure = nowCapPressure;
    }
    Vector3d & wallCapForce   = wallCap->getForceNonConst();
    Vector3d & wallCapMoment  = wallCap->getMomentNonConst();
    MPI_Allreduce(MPI_IN_PLACE, stressVoigt.data(), 6, MPI_DOUBLE, MPI_SUM, BOUNDARY_COMM);
    MPI_Allreduce(MPI_IN_PLACE, wallCapForce.data(), 3, MPI_DOUBLE, MPI_SUM, BOUNDARY_COMM);
    MPI_Allreduce(MPI_IN_PLACE, wallCapMoment.data(), 3, MPI_DOUBLE, MPI_SUM, BOUNDARY_COMM);
    if(rank == USED){
      stressVoigt = stressVoigt/wallBalls->findVolume()/scalingFactor;
      cout<<"current stress is: "<<stressVoigt(0)<<" "<<stressVoigt(1)<<" "<<stressVoigt(2)<<" "<<stressVoigt(3)<<" "<<stressVoigt(4)<<" "<<stressVoigt(5)<<endl;
    }

    /*########################################################################*/
    /*############# EXCHANGE FOR MEMBRANE INTERACTION (<-RECV) ###############*/
    /*########################################################################*/
    start_time4 = std::chrono::steady_clock::now();
    if(step % T == 0){
      if(rank < USED){
        MPI_Status r_status; MPI_Recv(&recv, 1, MPI_INT, USED, rank, BOUNDARY_COMM, &r_status);
        if(recv != 0){
          MPI_Status delivery_status; boundaryDataRecv.resize(recv);
          MPI_Recv(boundaryDataRecv.data(), recv, BOUNDARY_RECV, USED, rank, BOUNDARY_COMM, &delivery_status);
          for(int i = 0; i < recv; ++i){
            int gid = boundaryDataRecv[i]._id;
            grainsWorld[gid]->addForce(boundaryDataRecv[i]._force);
            grainsWorld[gid]->addMoment(boundaryDataRecv[i]._moment);
          }
          boundaryDataRecv.clear(); boundaryDataRecv.shrink_to_fit();
        }
      }else{
        for(int i = 0; i < USED; ++i){
          MPI_Request request; send = totalBoundaryDataRecv[i].size();
          MPI_Isend(&send, 1, MPI_INT, i, i, BOUNDARY_COMM, &request);
          if(send != 0){ MPI_Send(totalBoundaryDataRecv[i].data(), send, BOUNDARY_RECV, i, i, BOUNDARY_COMM); }
          totalBoundaryDataRecv[i].clear(); totalBoundaryDataRecv[i].shrink_to_fit();
          MPI_Wait(&request, MPI_STATUS_IGNORE);
        }
      }
    }
    end_time4 = std::chrono::steady_clock::now(); diff4 = end_time4 - start_time4; seconds4 = diff4.count();
    T_comm += seconds4;

    /*########################################################################*/
    /*############## MOVE WALLCAP AND WALLPLANE STRAIN CONTROL ###############*/
    /*########################################################################*/
    if(step < tConsolidate){
      wallCap->takeTimeStepConsolidate(pressure, dt);
    }else{
      wallCap->takeTimeStepStrainControl(Vector3d(0.,0.,dHeight),dt);
    }
    wallPlane->takeTimeStepStrainControl(Vector3d(0.,0.,0.),dt);

    /*########################################################################*/
    /*############### TAKE TIME STEP FOR GRAINS AND MEMBRANE #################*/
    /*########################################################################*/
    auto start_time3 = std::chrono::steady_clock::now();
    if(rank < USED){
      for(int i = 0; i < num_grains; ++i){
        if(belongToThisRank[i]==1){
          find_pos(grainsWorld[i]->getPosition(), _idx, _idy, _idz);
          _old_bin_id = to_bin_id(_idx, _idy, _idz);
          grainsWorld[i]->takeTimestep(gDamping, dt);
          find_pos(grainsWorld[i]->getPosition(), _idx, _idy, _idz);
          _proc_id = to_proc_id(_idx, _idy, _idz);
          _bin_id = to_bin_id(_idx, _idy, _idz);
          if(_proc_id == rank){
            if(_bin_id != _old_bin_id){
              /* DELETE */
              if(firstGrainIdList[_old_bin_id] == i){
                firstGrainIdList[_old_bin_id] = grainIdList[i];
              }else{
                j = firstGrainIdList[_old_bin_id];
                while(grainIdList[j] != i) j = grainIdList[j];
                grainIdList[j] = grainIdList[i];
              }
              grainIdList[i] = -1;
              bins[_old_bin_id].erase(i);
              /* ADD */
              if(firstGrainIdList[_bin_id] == -1){
                firstGrainIdList[_bin_id] = i;
              }else{
                j = firstGrainIdList[_bin_id];
                while(grainIdList[j] != -1) j = grainIdList[j];
                grainIdList[j] = i;
              }
              binIdList[i] = _bin_id;
              bins[_bin_id].insert(i);
            }
          }
          else{
            packer[_proc_id].emplace_back(i, grainsWorld[i]->getNodeShears().size(), grainsWorld[i]->getQuat(),
              grainsWorld[i]->getPosition(), grainsWorld[i]->getVelocity(), grainsWorld[i]->getOmega());
            packerNodeShears[_proc_id].insert(packerNodeShears[_proc_id].end(),
              grainsWorld[i]->getNodeShears().begin(),grainsWorld[i]->getNodeShears().end());
            packerNodeContact[_proc_id].insert(packerNodeContact[_proc_id].end(),
              grainsWorld[i]->getNodeContact().begin(),grainsWorld[i]->getNodeContact().end());
            packerNodeNormals[_proc_id].insert(packerNodeNormals[_proc_id].end(),
              grainsWorld[i]->getNodeNormals().begin(),grainsWorld[i]->getNodeNormals().end());
            /* DELETE */
            if(firstGrainIdList[_old_bin_id] == i){
              firstGrainIdList[_old_bin_id] = grainIdList[i];
            }else{
              j = firstGrainIdList[_old_bin_id];
              while(grainIdList[j] != i) {j = grainIdList[j]; }
              grainIdList[j] = grainIdList[i];
            }
            grainIdList[i] = -1;
            binIdList[i] = -1;
            bins[_old_bin_id].erase(i);
            belongToThisRank[i] = 0;
          }
        }
        else if(belongToThisRank[i] == 2) {
          grainIdList[i] = -1;
          binIdList[i] = -1;
          belongToThisRank[i] = 0;
        }
      }
    }else{
      //if(step % T == 0){ wallBalls->forceBetweenWalls(externalWall);}
      wallBalls->rotateTopLayer(wallCap);
      wallBalls->takeTimeStep(dt, "internalWall");
      externalWall->rotateTopLayer(wallCap);
      externalWall->takeTimeStep(dt, "externalWall");

      if(step >= tConsolidate) {
          wallBalls->moveHeight(dHeight);
          externalWall->moveHeight(dHeight);
      }
    }

    /* particle migration */
    if(USED > 1 && rank < USED){ grainsExchange(indir); }
    auto end_time3 = std::chrono::steady_clock::now(); std::chrono::duration<double> diff3 = end_time3 - start_time3; double seconds3 = diff3.count();
    T_migrate += seconds3;

    /*########################################################################*/
    /*################ OUTPUT SIMULATIOIN RESULT FOR DEBUG ###################*/
    /*########################################################################*/
    double kinematicEnergy = 0.0;
    double maxVelocity = 0.0;
    double maxOmega = 0.0;
    int    totalCN = 0;
    if(rank < USED){
      for(int i = 0; i < num_grains; ++i){
        if(belongToThisRank[i] == 1){
          maxVelocity = max(maxVelocity, grainsWorld[i]->getVelocity().norm());
          maxOmega    = max(maxOmega, grainsWorld[i]->getOmega().norm());
          totalCN    += grainsWorld[i]->getNContact();
          kinematicEnergy += grainsWorld[i]->computeKineticEnergy();
        }
      }
    }
    // MPI calls  sendbuf       recvbuff          count     type       op       comm
    MPI_Allreduce(MPI_IN_PLACE, &kinematicEnergy, 1, MPI_DOUBLE, MPI_SUM, BOUNDARY_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &totalCN, 1, MPI_INT, MPI_SUM, BOUNDARY_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &maxVelocity, 1, MPI_DOUBLE, MPI_MAX, BOUNDARY_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &maxOmega, 1, MPI_DOUBLE, MPI_MAX, BOUNDARY_COMM);
    double averageCN = (double)totalCN/num_grains;

    if(rank == 0){
      Vector3d capNormal = wallCap->getNormalGlobal();
      cout<<"Max velocity is: "<<maxVelocity<<" Max omega is: "<<maxOmega<<" Cap normal is: "<<capNormal(0)<<" "<<capNormal(1)<<" "<<capNormal(2)<<endl;
    }
    if(step % 10 == 0){
      if(rank == 0){
        cout<<"THIS IS ITERATION: "<<step<<" kinematic energy is: "<<kinematicEnergy*scalingFactor*scalingFactor<<" average coordination number: "<<averageCN<<" current dHeight: "<<dHeight<<endl;
        cout<<"Normal processor: T_force: "<<T_force<<" T_border: "<<T_border<<" T_migrate: "<<T_migrate<<" T_comm: "<<T_comm<<endl;
      }
      if(rank == USED){
        double grainVolume = totalMass/grainDensity;
        double totalVolume = wallBalls->findVolume();
        double voidVolume = totalVolume - grainVolume;
        double voidRatio = voidVolume/grainVolume;
        cout<<"Membrane processor: T_force: "<<T_force<<" T_border: "<<T_border<<" T_migrate: "<<T_migrate<<" T_comm: "<<T_comm<<endl;
        cout<<"Volume is: "<<wallBalls->findVolume()<<" voidRatio is: "<<voidRatio<<" stress_33 is: "<<stressVoigt(2)<<" stress_11 is: "<<stressVoigt(0)<<endl;
        cout<<"############## END "<<step<<" ITERATIONS ################"<<endl<<endl;
      }
    }

    /*########################################################################*/
    /*########################### OUTPUT TO FILES ############################*/
    /*########################################################################*/
    /* output specimen volume and average coordination number */
    if(rank == USED && step % outputFrequency == 0){
      fprintf(volumeCNFile, "%.8f %.4f\n", wallBalls->findVolume(), averageCN);
      fflush(volumeCNFile);
    }

    if(rank < USED && outputForceChain && step % outputFrequency == 0){
      /* allreduce chain numbers */
      int chainNum = forceChain.size();
      MPI_Allreduce(MPI_IN_PLACE, &chainNum, 1, MPI_INT, MPI_SUM, CART_COMM);
      if(rank > 0){
        MPI_Request request;
        int sendNumber = forceChain.size();
        MPI_Isend(&sendNumber,1,MPI_INT,0,rank,CART_COMM,&request);
        if(sendNumber != 0){
          MPI_Send(posChain.data(),sendNumber*6,MPI_DOUBLE,0,rank,CART_COMM);
          MPI_Send(forceChain.data(),sendNumber*3,MPI_DOUBLE,0,rank,CART_COMM);
        }
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }
      if(rank == 0){
        FILE * chainPositionFile  = fopen((outdir+"posChain_"+testName+"_"+std::to_string(step/outputFrequency)+".dat").c_str(), "w");
        FILE * chainForceFile  = fopen((outdir+"forceChain_"+testName+"_"+std::to_string(step/outputFrequency)+".dat").c_str(), "w");
        fprintf(chainPositionFile, "%d\n", chainNum);
        fprintf(chainForceFile, "%d\n", chainNum);
        for(int i = 0; i < forceChain.size(); ++i){
          Vector6d position = posChain[i];
          Vector3d force    = forceChain[i];
          fprintf(chainPositionFile, "%.4f %.4f %.4f %.4f %.4f %.4f\n", position(0), position(1), position(2), position(3), position(4), position(5));
          fprintf(chainForceFile, "%.4f %.4f %.4f\n", force(0), force(1), force(2));
        }
        for(int _rank = 1; _rank < USED; ++_rank){
          MPI_Status status;
          int recvNum = 0;
          MPI_Recv(&recvNum,1,MPI_INT,_rank,_rank,CART_COMM,&status);
          if(recvNum != 0){
            MPI_Status delivery_status;
            vector<Vector6d> positions;
            vector<Vector3d> forces;
            positions.resize(recvNum);
            forces.resize(recvNum);
            MPI_Recv(positions.data(),recvNum*6,MPI_DOUBLE,_rank,_rank,CART_COMM,&delivery_status);
            MPI_Recv(forces.data(),recvNum*3,MPI_DOUBLE,_rank,_rank,CART_COMM,&delivery_status);
            for(int i = 0; i < recvNum; ++i){
              Vector6d position = positions[i];
              Vector3d force    = forces[i];
              fprintf(chainPositionFile, "%.4f %.4f %.4f %.4f %.4f %.4f\n", position(0), position(1), position(2), position(3), position(4), position(5));
              fprintf(chainForceFile, "%.4f %.4f %.4f\n", force(0), force(1), force(2));
            }
          }
        }
        fclose(chainPositionFile);
        fclose(chainForceFile);
      }
    }
  }
}
#endif /*MPI_SOLVER_H_*/
