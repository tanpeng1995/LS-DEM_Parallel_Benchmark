#ifndef SOLVER_HPP_
#define SOLVER_HPP_
#include <algorithm>
#include "collision_resolution.hpp"
#include "Grain3d.hpp"
#include "WallPlane.hpp"
#include "WallBalls.hpp"
#include "WallTopography.hpp"
#include "init_simulation.hpp"

void simulate_one_step( double dt, int step, string indir){
  for(int i = 0; i < num_grains; ++i){
    /* Apply gravity force */
    /*##########################################################*/
    grainsWorld[i]->addVelocityAcceleration(Vector3d(0.,0.,-9.81)/scalingFactor,dt);
    Vector3d vel = grainsWorld[i]->getVelocity();
    Vector3d omg = grainsWorld[i]->getOmega();
    /* Zero out bad grains */
    if(isnan(vel.norm()) || isnan(omg.norm()) || isnan(-vel.norm()) || isnan(-omg.norm()) ){
      vel = Vector3d(0.,0.,0.);
      omg = Vector3d(0.,0.,0.);
    }
    /* Apply damping */
    grainsWorld[i]->changeVelocity( (1. - dt * gDamping ) * vel );
    grainsWorld[i]->changeOmega( (1. - dt * gDamping ) * omg );
    /*##########################################################*/
  }
  //cout<<"finish apply gravity force at step: "<<step<<endl;

  /* Force resolution */
  vector<vector<int>> grainContactId(num_grains,vector<int>());
  vector<vector<ContactInfo>> grainContactInfo(num_grains, vector<ContactInfo>());
  /* Membrane force, it changes velocity */
  for(int i = 0; i < num_grains; ++i){
    //wallBalls[0]->findWallForceMoment(grainsWorld[i], grainContactInfo[i]);
    //wallBalls[1]->findWallForceMoment(grainsWorld[i], grainContactInfo[i]);
  }
  /* inter-grain force */
  for(int i = 0; i < num_grains; ++i){
    for(int j = i+1; j < num_grains; ++j){
      if(grainsWorld[i]->bCircleCheck(grainsWorld[j])){
        grainsWorld[i]->findInterGrainContact(grainsWorld[j],grainContactInfo[i],grainContactId[j],grainContactId[i]);
      }
    }
    if(wallTopography->bCircleCheck(grainsWorld[i])){
      wallTopography->findGrainWallContactImpulse(grainsWorld[i], dt, grainContactInfo[i]);
    }
    if(wallPlanes[LEFT]->bCircleCheck(grainsWorld[i])){
      wallPlanes[LEFT]->findGrainWallContactImpulse(grainsWorld[i], dt, grainContactInfo[i]);
    }
    if(wallPlanes[RIGHT]->bCircleCheck(grainsWorld[i])){
      wallPlanes[RIGHT]->findGrainWallContactImpulse(grainsWorld[i], dt, grainContactInfo[i]);
    }
    if(wallPlanes[UP]->bCircleCheck(grainsWorld[i])){
      wallPlanes[UP]->findGrainWallContactImpulse(grainsWorld[i], dt, grainContactInfo[i]);
    }
    if(wallPlanes[DOWN]->bCircleCheck(grainsWorld[i])){
      wallPlanes[DOWN]->findGrainWallContactImpulse(grainsWorld[i], dt, grainContactInfo[i]);
    }
    if(grainsWorld[i]->getPosition()(0) <= 1000. && wallPlanes[FRONT]->bCircleCheck(grainsWorld[i])){
      wallPlanes[FRONT]->findGrainWallContactImpulse(grainsWorld[i], dt, grainContactInfo[i]);
    }
    if(grainsWorld[i]->getPosition()(0) <= 1000. && wallPlanes[BACK]->bCircleCheck(grainsWorld[i])){
      wallPlanes[BACK]->findGrainWallContactImpulse(grainsWorld[i], dt, grainContactInfo[i]);
    }
  }
  //cout<<"finish inter-grain force at step: "<<step<<endl;

  /* Graph partition */
  std::vector<int> degreeList(num_grains,0);
  std::vector<bool> checkList(num_grains,false);
  std::priority_queue<DegreeQueue, std::vector<DegreeQueue>, CompareDegree> iterateQueue;
  int islandsCount = 0;
  std::vector<std::set<int, CompareSet>> islands;
  for(int i = 0; i < num_grains; ++i){
    degreeList[i] = grainContactId[i].size();
    iterateQueue.emplace(i, grainContactId[i].size());
  }
  /* The CutHill-McKee Algorithm */
  //initialize sets for contact island and priority_queue for search
  int totalCheckedGrainNumber = 0;
  std::priority_queue<DegreeQueue, std::vector<DegreeQueue>, CompareDegree> degreeQueue;
  while(totalCheckedGrainNumber < num_grains){
    //find the grain with minimum degree number
    int minDegree = std::numeric_limits<int>::max();
    int minDegreeId = 0;
    DegreeQueue firstPop = iterateQueue.top();
    while(checkList[firstPop._id]){
      iterateQueue.pop();
      firstPop = iterateQueue.top();
    }
    minDegree = firstPop._degree;
    minDegreeId = firstPop._id;
    iterateQueue.pop();

    DegreeQueue firstDegree;
    //add minDegreeId into set, mark as checked and add its neighbors into degreeQueue
    islands.push_back(std::set<int, CompareSet>{minDegreeId});
    checkList[minDegreeId] = true;
    for(int i = 0; i < grainContactId[minDegreeId].size(); ++i){
      int gid = grainContactId[minDegreeId][i];
      int degree = degreeList[gid];
      degreeQueue.emplace(gid, degree);
    }
    //while degreeQueue is not empty
    while(!degreeQueue.empty()){
      //pop the first item, if is not in set, then add in set
      firstDegree = degreeQueue.top();
      degreeQueue.pop();
      int minId = firstDegree._id;
      auto got = islands[islandsCount].find(minId);
      if( got == islands[islandsCount].end() ){ // if not in set
        islands[islandsCount].insert(minId);    // insert
        checkList[minId] = true;                // mark true
        for(int i = 0; i < grainContactId[minId].size(); ++i){
          int gid = grainContactId[minId][i];
          int degree = degreeList[gid];
          degreeQueue.emplace(gid, degree);
        }
      }
    } // finish one island
    totalCheckedGrainNumber += islands[islandsCount].size(); // count totalnumber
    islandsCount ++; // islandsCount increase by 1
  }
  cout<<"number of islands: "<<islandsCount<<endl;

  std::vector<vector<int>> vectorizedSet(islandsCount);
  for(int i = 0; i < islandsCount; ++i){
    vectorizedSet[i].assign(islands[i].begin(), islands[i].end());
  }

  /* build contact lists */
  std::vector<vector<ContactInfo>> contactLists;
  for(int i = 0; i < islandsCount; ++i){
    std::vector<ContactInfo> cList;
    for(int j = 0; j < vectorizedSet[i].size(); ++j){
      int gid = vectorizedSet[i][j];
      cList.insert(cList.end(), grainContactInfo[gid].begin(), grainContactInfo[gid].end());
    }
    contactLists.push_back(cList);
  }
  //cout<<"finish contact island at step: "<<step<<endl;

  /* iteratively solve */
  const int tempN = 3;
  for(int tempStep = 0; tempStep < tempN; tempStep++){
    /* collision resolution */
    for(int i = 0; i < contactLists.size(); ++i){
      collision_resolution(contactLists[i], step);
    }
    /* take temp time step */
    for(int i = 0; i < num_grains; ++i){
      grainsWorld[i]->takeTempTimestep(dt);
    }
    /* update grains positions, correct velocities */
    for(int i = 0; i < num_grains; ++i){
      if(wallPlanes[LEFT]->bTempCircleCheck(grainsWorld[i])){
        wallPlanes[LEFT]->findGrainWallContact(grainsWorld[i], dt);
      }
      if(wallPlanes[RIGHT]->bTempCircleCheck(grainsWorld[i])){
        wallPlanes[RIGHT]->findGrainWallContact(grainsWorld[i], dt);
      }
      if(wallPlanes[UP]->bTempCircleCheck(grainsWorld[i])){
        wallPlanes[UP]->findGrainWallContact(grainsWorld[i], dt);
      }
      if(wallPlanes[DOWN]->bTempCircleCheck(grainsWorld[i])){
        wallPlanes[DOWN]->findGrainWallContact(grainsWorld[i], dt);
      }
      if(grainsWorld[i]->getTempPosition()(0) <= 1000. && wallPlanes[FRONT]->bTempCircleCheck(grainsWorld[i])){
        wallPlanes[FRONT]->findGrainWallContact(grainsWorld[i], dt);
      }
      if(grainsWorld[i]->getTempPosition()(0) <= 1000. && wallPlanes[BACK]->bTempCircleCheck(grainsWorld[i])){
        wallPlanes[BACK]->findGrainWallContact(grainsWorld[i], dt);
      }
      if(wallTopography->bTempCircleCheck(grainsWorld[i])){
        wallTopography->findGrainWallContact(grainsWorld[i], dt);
      }
    }
    /* reset relative velocity */
    for(int i = 0; i < contactLists.size(); ++i){
      for(auto & item : contactLists[i]){
        if(item._slave < 0){
          item.computeRelVelocityWithTopographyThreshold();
        }else if(item._slave >= num_grains){
          const int wallId = whichWallBall(item._slave);
          item.computeRelVelocityWithBallsThreshold(wallBalls[wallId]->getVelocity(item._slave - wallBalls[wallId]->getStartIdx()));
        }else{
          item.computeRelVelocityThreshold();
        }
      }
    }
  }
  /* Take time step */
  //wallBalls[0]->takeTimeStep(dt);
  //wallBalls[0]->correctMembraneVelocity(dt);
  //wallBalls[1]->takeTimeStep(dt);
  //wallBalls[1]->correctMembraneVelocity(dt);

  for(int i = 0; i < num_grains; ++i){
    Vector3d vel = grainsWorld[i]->getVelocity();
    if(vel.norm() > 400.){
      grainsWorld[i]->changeVelocity(vel/vel.norm()*300.);
    }
    Vector3d omg = grainsWorld[i]->getOmega();
    if(omg.norm() > 50.){
      grainsWorld[i]->changeOmega(omg/omg.norm()*50.);
    }
    grainsWorld[i]->takeTimestep(dt);
  }
  double kinematicEnergy = 0.0;
  double maxVelocity = 0.0;
  double maxOmg = 0.0;
  for(int i = 0; i < num_grains; ++i){
    kinematicEnergy += grainsWorld[i]->computeKineticEnergy();
    maxVelocity = max( maxVelocity, grainsWorld[i]->getVelocity().norm() );
    maxOmg = max( maxOmg, grainsWorld[i]->getOmega().norm() );
  }
  if(step % 10 == 0){ cout<<"THIS IS ITERATION: "<<step<<" kinematic energy is: "<<kinematicEnergy<<" max velocity is: "<<maxVelocity<<" max Omega is: "<<maxOmg<<endl; }
}
#endif /*SOLVER_HPP_*/
