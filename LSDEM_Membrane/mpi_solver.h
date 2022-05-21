/*
 * MPI_SOLVER_H_
 * Created on: June 21, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_SOLVER_H_
#define MPI_SOLVER_H_
#include "mpi_border_communication.h"
#include "mpi_particle_migration.h"

const int T = 5;

void simulate_one_step(int rank, double dt, int step, string indir, string outdir, FILE * volumeCNFile, double & T_force, double & T_border, double & T_migrate, double & T_comm){
  if(rank < USED){
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
    if(USED > 1){ updateRightLeftSide(indir); updateBackFrontSide(indir); updateDownUpSide(indir); }
    clearForces(); updateAllPoints();
    auto end_time1 = std::chrono::steady_clock::now(); std::chrono::duration<double> diff1 = end_time1 - start_time1; double seconds1 = diff1.count();
    T_border += seconds1;

    /*########################################################################*/
    /*########################## FORCE RESOLUTION ############################*/
    /*########################################################################*/
    vector<Vector6d> posChain;
    vector<Vector3d> forceChain;
    auto start_time2 = std::chrono::steady_clock::now();
    int j, _idx, _idy, _idz, _bin_id, _proc_id, _old_bin_id;
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
        /*###################### CONTACT WITH MEMBRANE #####################*/
        if(step % T == 0){
          //if(wallBalls->bCircleCheck(grainsWorld[i])){
            wallBalls->findWallForceMoment(grainsWorld[i],dt,stressVoigt);
          //}
          //if(externalWall->bCircleCheck(grainsWorld[i])){
            externalWall->findWallForceMoment(grainsWorld[i],dt,stressVoigt);
          //}
        }
        /* end boundary interactions */
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
      /* end force resolution */
    }
    auto end_time2 = std::chrono::steady_clock::now(); std::chrono::duration<double> diff2 = end_time2 - start_time2; double seconds2 = diff2.count();
    T_force += seconds2;

    /*########################################################################*/
    /*############### FIND CAP AND PLANE FORCES WITH GRAINS ##################*/
    /*########################################################################*/
    Vector3d & wallCapForce   = wallCap->getForceNonConst();
    Vector3d & wallCapMoment  = wallCap->getMomentNonConst();
    Vector3d & wallPlaneForce = wallPlane->getForceNonConst();
    MPI_Allreduce(MPI_IN_PLACE, wallCapForce.data(), 3, MPI_DOUBLE, MPI_SUM, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, wallCapMoment.data(), 3, MPI_DOUBLE, MPI_SUM, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, wallPlaneForce.data(), 3, MPI_DOUBLE, MPI_SUM, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, stressVoigt.data(), 6, MPI_DOUBLE, MPI_SUM, CART_COMM);
    nowCapPressure = abs( wallCapForce.dot(wallCap->getNormalGlobal()) )/ wallCap->getArea();
    planePressure  = abs( wallPlaneForce.dot(wallPlane->getNormalGlobal()) )/ wallCap->getArea();
    stressVoigt = stressVoigt/wallBalls->findVolume()/scalingFactor;
    if(rank == 0){
      cout<<"current cap pressure is: "<<nowCapPressure/scalingFactor<<" last cap pressure is: "<<prevCapPressure/scalingFactor<<" planePressure is: "<<planePressure/scalingFactor<<endl;
      cout<<"current stress is: "<<stressVoigt(0)<<" "<<stressVoigt(1)<<" "<<stressVoigt(2)<<" "<<stressVoigt(3)<<" "<<stressVoigt(4)<<" "<<stressVoigt(5)<<endl;
    }
    prevCapPressure = nowCapPressure;

    /*########################################################################*/
    /*############### TAKE TIME STEP FOR GRAINS AND MEMBRANE #################*/
    /*########################################################################*/
    auto start_time3 = std::chrono::steady_clock::now();
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
    /*########################################################################*/
    /*################ BROADCAST BALLWALL POSITION VELOCITY ##################*/
    /*########################################################################*/
    /* inner wallBalls */
    vector<int> ballsRank = wallBalls->getBallsRank();
    vector<Vector3d> ballsFor = wallBalls->getForces();
    vector<MPIBalls> mpiBalls;
    for(int i = 0; i < ballsRank.size(); ++i){
      int ball_id = ballsRank[i];
      mpiBalls.emplace_back(rank, ball_id, ballsFor[ball_id]);
    }
    int nBalls = mpiBalls.size();
    vector<int> buffer(USED,0);
    MPI_Allgather(&nBalls, 1, MPI_INT, buffer.data(), 1, MPI_INT, CART_COMM);
    int allBalls = 0;
    for(int i = 0; i < USED; ++i){ allBalls += buffer[i]; }
    vector<int> displacement(USED,0);
    for(int i = 1; i < USED; ++i){ displacement[i] = displacement[i-1] + buffer[i-1]; }
    vector<MPIBalls> mpiAllBalls(allBalls);
    MPI_Allgatherv(mpiBalls.data(), buffer[rank], MPI_BALL, mpiAllBalls.data(), buffer.data(), displacement.data(), MPI_BALL, CART_COMM);
    for(auto & item : mpiAllBalls){
      wallBalls->changeForce(item._id, item._force);
    }

    /* output wallBalls in each island for visualization */
    if(outputBallRank && step % outputFrequency == 0 && rank == 0){
      FILE * ballRankFile  = fopen((outdir+"ballRank_"+testName+"_"+std::to_string(step/outputFrequency)+".dat").c_str(), "w");
      fprintf(ballRankFile, "%d\n", mpiAllBalls.size());
      vector<MPIBalls> tempMPIAllBalls;
      for(int i = 0; i < USED-1; ++i){
        tempMPIAllBalls.assign(mpiAllBalls.begin()+displacement[i], mpiAllBalls.begin()+displacement[i+1]);
        for(auto & item : tempMPIAllBalls){
          fprintf(ballRankFile, "%d %d\n", item._id, i);
        }
      }
      tempMPIAllBalls.assign(mpiAllBalls.begin()+displacement[USED-1], mpiAllBalls.end());
      for(auto & item : tempMPIAllBalls){
        fprintf(ballRankFile, "%d %d\n", item._id, USED-1);
      }
      fflush(ballRankFile);
    }

    /*external wallBalls*/

    ballsRank = externalWall->getBallsRank();
    ballsFor = externalWall->getForces();
    mpiBalls.clear(); mpiBalls.shrink_to_fit();
    mpiAllBalls.clear(); mpiAllBalls.shrink_to_fit();
    for(int i = 0; i < ballsRank.size(); ++i){
      int ball_id = ballsRank[i];
      mpiBalls.emplace_back(rank, ball_id, ballsFor[ball_id]);
    }
    nBalls = mpiBalls.size();
    MPI_Allgather(&nBalls, 1, MPI_INT, buffer.data(), 1, MPI_INT, CART_COMM);
    allBalls = 0;
    for(int i = 0; i < USED; ++i){ allBalls += buffer[i]; }
    for(int i = 0; i < USED; ++i){ displacement[i] = 0; }
    for(int i = 1; i < USED; ++i){ displacement[i] = displacement[i-1] + buffer[i-1]; }
    mpiAllBalls.resize(allBalls);
    MPI_Allgatherv(mpiBalls.data(), buffer[rank], MPI_BALL, mpiAllBalls.data(), buffer.data(), displacement.data(), MPI_BALL, CART_COMM);
    for(auto & item : mpiAllBalls){
      externalWall->changeForce(item._id, item._force);
    }

    /*########################################################################*/
    /*#################### MOVE ALL WALLS STRAIN CONTROL #####################*/
    /*########################################################################*/
    
      //wallBalls->computeInternalForceOnTopLayer(wallCap, dt); // may cause instability
      //externalWall->computeInternalForceOnTopLayer(wallCap, dt); // may cause instability
      if(step < tConsolidate){
        wallCap->takeTimeStepConsolidate(pressure, dt);
      }else{
        wallCap->takeTimeStepStrainControl(Vector3d(0.,0.,dHeight),dt);
      }
      wallPlane->takeTimeStepStrainControl(Vector3d(0.,0.,0.),dt);
      if(step > tConsolidate){
        wallBalls->moveHeight(dHeight);
        externalWall->moveHeight(dHeight);
      }
      wallBalls->rotateTopLayer(wallCap);
      wallBalls->takeTimeStep(dt, "internalWall");
      externalWall->rotateTopLayer(wallCap);
      externalWall->takeTimeStep(dt, "externalWall");

    /*########################################################################*/
    /*################## CROSS-BLOCK MIGRATION FOR GRAIN #####################*/
    /*########################################################################*/
    if(USED > 1){ grainsExchange(indir); }
    auto end_time3 = std::chrono::steady_clock::now(); std::chrono::duration<double> diff3 = end_time3 - start_time3; double seconds3 = diff3.count();
    T_migrate += seconds3;

    /*########################################################################*/
    /*################ OUTPUT SIMULATIOIN RESULT FOR DEBUG ###################*/
    /*########################################################################*/
    double kinematicEnergy = 0.0;
    double maxVelocity = 0.0;
    double maxOmega = 0.0;
    int    totalCN = 0;
    for(int i = 0; i < num_grains; ++i){
      if(belongToThisRank[i] == 1){
        maxVelocity = max(maxVelocity, grainsWorld[i]->getVelocity().norm());
        maxOmega    = max(maxOmega, grainsWorld[i]->getOmega().norm());
        totalCN    += grainsWorld[i]->getNContact();
        kinematicEnergy += grainsWorld[i]->computeKineticEnergy();
      }
    }
    // MPI calls  sendbuf       recvbuff          count     type       op       comm
    MPI_Allreduce(MPI_IN_PLACE, &kinematicEnergy, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &totalCN, 1, MPI_INT, MPI_SUM, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &maxVelocity, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &maxOmega, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    double averageCN = (double)totalCN/num_grains;

    if(rank == 0){
      Vector3d capNormal = wallCap->getNormalGlobal();
      cout<<"Max velocity is: "<<maxVelocity<<" Max omega is: "<<maxOmega<<" Cap normal is: "<<capNormal(0)<<" "<<capNormal(1)<<" "<<capNormal(2)<<endl;
    }

    if(step % 10 == 0){
      if(rank == 0){
        double grainVolume = totalMass/grainDensity;
        double totalVolume = wallBalls->findVolume();
        double voidVolume = totalVolume - grainVolume;
        double voidRatio = voidVolume/grainVolume;
        cout<<"THIS IS ITERATION: "<<step<<" kinematic energy is: "<<kinematicEnergy*scalingFactor*scalingFactor<<" average coordination number: "<<averageCN<<" current dHeight: "<<dHeight<<endl;
        cout<<"Normal processor: T_force: "<<T_force<<" T_border: "<<T_border<<" T_migrate: "<<T_migrate<<" T_comm: "<<T_comm<<endl;
        cout<<"Volume is: "<<wallBalls->findVolume()<<" voidRatio is: "<<voidRatio<<" stress_33 is: "<<stressVoigt(2)<<" stress_11 is: "<<stressVoigt(0)<<endl;
        cout<<"############## END "<<step<<" ITERATIONS ################"<<endl<<endl;
      }
    }

    /*########################################################################*/
    /*########################### OUTPUT TO FILES ############################*/
    /*########################################################################*/
    /* output specimen volume and average coordination number */
    if(step % outputFrequency == 0){
      fprintf(volumeCNFile, "%.8f %.4f\n", wallBalls->findVolume(), averageCN);
      fflush(volumeCNFile);
    }

    if(outputForceChain && step % outputFrequency == 0){
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
    /* end output */
  }
}
#endif /*MPI_SOLVER_H_*/
