/*
 * MPI_CYLINDER_WALL_SOLVER_H_
 * Created on: June 21, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_CYLINDER_WALL_SOLVER_H_
#define MPI_CYLINDER_WALL_SOLVER_H_
#include "mpi_border_communication.h"
#include "mpi_particle_migration.h"

const int T = 1;
void simulate_one_step(
  int rank, double dt, int step, string indir,
  double & T_force, double & T_border, double & T_migrate, double & T_update){
  /* set random seed */
  std::mt19937 gen(37*step+17);
  std::uniform_real_distribution<double> rand_real(-1.0, 1.0);
  stressVoigt << 0.,0.,0.,0.,0.,0.;

  if(rank < USED){
    //T_border
    MPI_Barrier(CART_COMM);
    auto start_time1 = std::chrono::steady_clock::now();
    if(USED > 1){ updateRightLeftSide(indir); updateBackFrontSide(indir); updateDownUpSide(indir);}
    auto end_time1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff1 = end_time1 - start_time1;
    double seconds1 = diff1.count();

    MPI_Barrier(CART_COMM);
    auto start_time4 = std::chrono::steady_clock::now();
    clearForces(step); updateAllPoints();
    auto end_time4 = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff4 = end_time4 - start_time4;
    double seconds4 = diff4.count();

    //T_force
    MPI_Barrier(CART_COMM);
    auto start_time2 = std::chrono::steady_clock::now();
    int j, _idx, _idy, _idz, _bin_id, _proc_id, _old_bin_id;
    #pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads) private(j, _idx, _idy, _idz, _bin_id, _proc_id, _old_bin_id) shared(grainsWorld, belongToThisRank, grainIdList, topPlane, bottomPlane, firstGrainIdList)
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
          binNum+binsInDimX*binsInDimY,              binNum+binsInDimX*binsInDimY+1,
          binNum+binsInDimX*binsInDimY+binsInDimX-1, binNum+binsInDimX*binsInDimY+binsInDimX,
          binNum+binsInDimX*binsInDimY+binsInDimX+1};

      if(belongToThisRank[i]==1){
        j = grainIdList[i];
        //inside the same bins
        while(j != -1){
          if(grainsWorld[i]->bCircleCheck(grainsWorld[j])){
            if(i < j){ grainsWorld[i]->findInterparticleForceMoment(grainsWorld[j],dt,stressVoigt); }
            else{ grainsWorld[j]->findInterparticleForceMoment(grainsWorld[i],dt,stressVoigt); }
          }
          j = grainIdList[j];
        }
        //consider 3, 6, 7, 8 in the same level as the bins of interest and 0~9 at upper level
        for (int k = 0; k < 13; ++k) {
          j = firstGrainIdList[binArray[k]];
          while(j != -1){
            if(grainsWorld[i]->bCircleCheck(grainsWorld[j])){
              if(i < j){ grainsWorld[i]->findInterparticleForceMoment(grainsWorld[j],dt,stressVoigt); }
              else{ grainsWorld[j]->findInterparticleForceMoment(grainsWorld[i],dt,stressVoigt); }
            }
            j = grainIdList[j];
          }
        }
        //grainsWorld[i]->applyAcceleration(Vector3d(rand_real(gen)*5.,rand_real(gen)*5.,rand_real(gen)*5.));
        /* bottom and top plane */
        if(bottomPlane->bCircleCheckNonPenetration(grainsWorld[i])){
          bottomPlane->findWallForceMoment(grainsWorld[i],dt);
        }
        if(topPlane->bCircleCheckNonPenetration(grainsWorld[i])){
          topPlane->findWallForceMoment(grainsWorld[i],dt);
        }
        /* wallCylinder */
        if(step % T == 0){
          if(wallCylinder->bCircleCheck(grainsWorld[i])){
            wallCylinder->findWallForceMomentNonFriction(grainsWorld[i],dt);
          }
        }

        /* wallBowl */
        /*
        if(step % T == 0){
          if(wallBowl->bCircleCheck(grainsWorld[i])){
            wallBowl->findWallForceMomentNonFriction(grainsWorld[i],dt);
          }
        }
        */

      }else if(belongToThisRank[i] == 2) {
        // 1. bin_id for ghosts 2. condition branch if target belongs to this rank
        // no calculations for inside the same bins
        for (int k = 0; k < 13; ++k) {
          find_pos(binArray[k], _idx, _idy, _idz);
          if(_idz < binsInDimZ-1 && _idz > 0 && _idy < binsInDimY-1 && _idy > 0 && _idx < binsInDimX-1 && _idx > 0) {
            j = firstGrainIdList[binArray[k]];
            while(j != -1){
              if(belongToThisRank[j] == 1 && grainsWorld[i]->bCircleCheck(grainsWorld[j])){
                if(i < j){ grainsWorld[i]->findInterparticleForceMoment(grainsWorld[j],dt,stressVoigt); }
                else{ grainsWorld[j]->findInterparticleForceMoment(grainsWorld[i],dt,stressVoigt); }
              }
              j = grainIdList[j];
            }
          }
        }
      }
    }

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
            //DELETE
            if(firstGrainIdList[_old_bin_id] == i){
              firstGrainIdList[_old_bin_id] = grainIdList[i];
            }else{
              j = firstGrainIdList[_old_bin_id];
              while(grainIdList[j] != i) j = grainIdList[j];
              grainIdList[j] = grainIdList[i];
            }
            grainIdList[i] = -1;
            bins[_old_bin_id].erase(i);
            //ADD
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
        }else{
          //cout<<"oh... particle migration idx: "<<i<<endl;
          packer[_proc_id].emplace_back(i, grainsWorld[i]->getNodeShears().size(), grainsWorld[i]->getQuat(), grainsWorld[i]->getPosition(), grainsWorld[i]->getVelocity(), grainsWorld[i]->getOmega());
          packerNodeShears[_proc_id].insert(packerNodeShears[_proc_id].end(), grainsWorld[i]->getNodeShears().begin(), grainsWorld[i]->getNodeShears().end());
          packerNodeContact[_proc_id].insert(packerNodeContact[_proc_id].end(), grainsWorld[i]->getNodeContact().begin(), grainsWorld[i]->getNodeContact().end());
          packerNodeNormals[_proc_id].insert(packerNodeNormals[_proc_id].end(), grainsWorld[i]->getNodeNormals().begin(), grainsWorld[i]->getNodeNormals().end());
          //DELETE
          if(firstGrainIdList[_old_bin_id] == i){
            firstGrainIdList[_old_bin_id] = grainIdList[i];
          }else{
            j = firstGrainIdList[_old_bin_id];
            while(grainIdList[j] != i) j = grainIdList[j];
            grainIdList[j] = grainIdList[i];
          }
          grainIdList[i] = -1;
          binIdList[i] = -1;
          bins[_old_bin_id].erase(i);
          belongToThisRank[i] = 0;
        }
      } else if(belongToThisRank[i] == 2) {
        //clear out border particles, since border_communication happens at every step
        grainIdList[i] = -1;
        binIdList[i] = -1;
        belongToThisRank[i] = 0;
      }
      /* end update grain*/
    }
    auto end_time2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff2 = end_time2 - start_time2;
    double seconds2 = diff2.count();

    //exchange grains
    MPI_Barrier(CART_COMM);
    auto start_time3 = std::chrono::steady_clock::now();
    if(USED > 1) { grainsExchange(indir); }
    auto end_time3 = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff3 = end_time3 - start_time3;
    double seconds3 = diff3.count();

    MPI_Allreduce(MPI_IN_PLACE, &seconds1, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &seconds2, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &seconds3, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, &seconds4, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    MPI_Allreduce(MPI_IN_PLACE, stressVoigt.data(), 6, MPI_DOUBLE, MPI_SUM, CART_COMM);

    T_border  += seconds1;
    T_force   += seconds2;
    T_migrate += seconds3;
    T_update  += seconds4;

    /*find highest pos*/
    double highPos = 0.;
    for(int i = 0; i < num_grains; ++i){
      if(belongToThisRank[i] == 1){
        vector<Vector3d> nodePositions = grainsWorld[i]->getPointList();
        for(auto pos : nodePositions){
          highPos = max(highPos, pos(2));
        }
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &highPos, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
    if(step <= 1000){
      //wallBowl->changeD(-0.025);
      //wallBowl->changeRadius(0.015);
      //wallBowl->changeHeight(0.02);
      //topPlane->increaseHeight(0.02);
      //wallCylinder->changeRadius(0.02);
    }
    double grainVolume = totalMass/grainDensity;
    double cylinderVolume = wallCylinder->getRadius()*wallCylinder->getRadius()*M_PI*highPos;
    double voidVolume = cylinderVolume - grainVolume;
    double voidRatio = voidVolume/grainVolume;

    /*print*/
    if(rank == 0){
      cout<<"At step "<<step<<" Maximum Grain Height: "<<highPos<<" topPlane heigh: "<<topPlane->getHeight()<<" wallCylinder radius: "<<wallCylinder->getRadius()<<" Void Ratio: "<<voidRatio<<endl;
    }
    if(rank == 0){
      stressVoigt = stressVoigt/cylinderVolume/scalingFactor;
      cout<<"current stress is: "<<stressVoigt(0)<<" "<<stressVoigt(1)<<" "<<stressVoigt(2)<<" "<<stressVoigt(3)<<" "<<stressVoigt(4)<<" "<<stressVoigt(5)<<endl;
    }

    if(step % 5 == 0){
      if(rank == 0){
        cout<<endl;
        cout<<" Force Interaction Time = " << T_force << " seconds. "<<" Border Communication Time = " << T_border << " seconds\n";
        cout<<" Particle Migration Time = " << T_migrate << " seconds. "<<" Info Update Time = " << T_update << " seconds\n";
        cout<<"############## END "<<step<<" ITERATIONS ################"<<endl<<endl;
      }
      /* end output */
    }
  }
}
#endif /*MPI_CYLINDER_WALL_SOLVER_H_*/
