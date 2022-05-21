#ifndef MPI_SOLVER_RIGID_HPP_
#define MPI_SOLVER_RIGID_HPP_
#include "mpi_border_communication.hpp"
#include "mpi_collision_resolution.hpp"
#include "mpi_move_grain.hpp"
#include "Grain3d.hpp"
#include "WallPlane.hpp"
#include "mpi_contact_normal.hpp"
#include "mpi_quick_hull.hpp"
#include "mpi_triangulation_pressure.hpp"
#include "mpi_delaunay_triangulation.hpp"
#include "algorithm"
#define START_LOOP(rank) if(rank < USED){
#define END_LOOP }
#define START_PARALLEL if(USED > 1){
#define END_PARALLEL }
#define BEGIN_MASTER_RANK(rank) if(rank == 0){
#define ELSE_MASTER_RANK else{
#define END_MASTER_RANK }
#define BEGIN_SLAVE_RANK(rank) if(rank != 0){
#define END_SLAVE_RANK }

void simulate_one_step(
  int rank, double dt, int step, string indir,
  double & T1, double & T2, double & T3, double & T4, double & T5, double & T6, double & T7, double & T8, double & T9 ){
  double t1 = 0., t2 = 0., t3 = 0., t4 = 0., t5 = 0., t6 = 0., t7 =0., t8 = 0., t9 = 0.;

  std::mt19937 gen(step*31+17);
  std::uniform_real_distribution<double> rand_real(-1.0, 1.0);
  if(rank < USED){
  for(int i = 0; i < num_grains; ++i){
    //if(step >= 100) grainsWorld[i].addVelocityAcceleration(Vector3d(0.,0.,-25.)/scalingFactor,dt);
    grainsWorld[i].addVelocityAcceleration(Vector3d(rand_real(gen),rand_real(gen),rand_real(gen))*5./scalingFactor,dt);
    Vector3d vel = grainsWorld[i].getVelocity();
    Vector3d omg = grainsWorld[i].getOmega();
    if(isnan(vel.norm()) || isnan(omg.norm()) || isnan(-vel.norm()) || isnan(-omg.norm()) ){
      vel << 0.,0.,0.;
      omg << 0.,0.,0.;
    }
    grainsWorld[i].changeVelocity((1. - dt * gDamping) * vel);
    grainsWorld[i].changeOmega((1. - dt * gDamping) * omg);
  }
  /*##########################################################*/
  /*################ BORDER COMMUNICATION ####################*/
  /*##########################################################*/
  MPI_Barrier(CART_COMM);
  auto start_time1 = std::chrono::steady_clock::now();
  if(USED > 1){
  updateRightLeftSide(indir);
  updateBackFrontSide(indir);
  updateDownUpSide(indir);
  }
  auto end_time1 = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff1 = end_time1 - start_time1;
  double seconds1 = diff1.count();
  t1 += seconds1;
  /*##########################################################*/
  /*################## FORCE RESOLUTION ######################*/
  /*##########################################################*/
  MPI_Barrier(CART_COMM);
  auto start_time2 = std::chrono::steady_clock::now();
  int j, _idx, _idy, _idz, _bin_id, _proc_id, _old_bin_id;
  /*grainContactId, which grains are contacting*/
  /*contact info*/
  vector<vector<int>> grainContactId(num_grains,vector<int>());
  vector<vector<ContactInfo>> grainContactInfo(num_grains, vector<ContactInfo>());
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
    vector<ContactInfo> singleContactList;
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
        if(grainsWorld[i].bCircleCheck(grainsWorld[j])){
          grainsWorld[i].findInterGrainContact(grainsWorld[j],grainContactInfo[i],grainContactId[i]);
        }
        j = grainIdList[j];
      }
      //consider 3, 6, 7, 8 in the same level as the bins of interest and 0~9 at upper level
      for (int k = 0; k < 13; ++k) {
        j = firstGrainIdList[binArray[k]];
        while(j != -1){
          if(grainsWorld[i].bCircleCheck(grainsWorld[j])){
            grainsWorld[i].findInterGrainContact(grainsWorld[j],grainContactInfo[i],grainContactId[i]);
          }
          j = grainIdList[j];
        }
      }
    }else if(belongToThisRank[i] == 2) {
      for (int k = 0; k < 13; ++k) {
        find_pos(binArray[k], _idx, _idy, _idz);
        if(_idz < binsInDimZ-1 && _idz > 0 && _idy < binsInDimY-1 && _idy > 0 && _idx < binsInDimX-1 && _idx > 0) {
          j = firstGrainIdList[binArray[k]];
          while(j != -1){
            if(belongToThisRank[j] == 1 && grainsWorld[i].bCircleCheck(grainsWorld[j])){
              grainsWorld[j].findInterGrainContact(grainsWorld[i],grainContactInfo[j],grainContactId[j]);
            }
            j = grainIdList[j];
          }
        }
      }
    }

  }
  auto end_time2 = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff2 = end_time2 - start_time2;
  double seconds2 = diff2.count();
  t2 += seconds2;
  /*##########################################################*/
  /*#################### GRAPH PARTITION #####################*/
  /*##########################################################*/
  //graph_partition
  MPI_Barrier(CART_COMM);
  auto start_time3 = std::chrono::steady_clock::now();
  int islandsCount = 0;
  std::vector<int> originRankList(num_grains,0);
  std::vector<int> coordinateNumberList(num_grains,0);
  std::vector<std::set<int, CompareSet>> islands;
  {
    if(rank != 0){
      /*#################### FOR NORMAL GRAINS #####################*/
      std::vector<int> grainIdsOfRank;              //list of grain of rank
      std::vector<int> contactIdCountOfEachGrain;   //list of #of contact for each grain
      std::vector<int> contactInfoCountOfEachGrain; //list of #of contact for each grain
      std::vector<int> allContactPairOfRank;        //list of all contacts of rank
      for(int i = 0; i < num_grains; ++i){
        if(belongToThisRank[i] == 1){
          grainIdsOfRank.push_back(i);
          contactIdCountOfEachGrain.push_back(grainContactId[i].size());
          contactInfoCountOfEachGrain.push_back(grainContactInfo[i].size());
          allContactPairOfRank.insert(allContactPairOfRank.end(), grainContactId[i].begin(), grainContactId[i].end());
        }
      }
      MPI_Request request;
      /*send number of grain in this rank and send number of contact in each grain*/
      int numOfGrainsInRank = grainIdsOfRank.size();
      MPI_Isend(&numOfGrainsInRank,1,MPI_INT,0,rank,CART_COMM,&request);
      if(numOfGrainsInRank != 0){
        MPI_Send(grainIdsOfRank.data(),numOfGrainsInRank,MPI_INT,0,rank,CART_COMM);
        MPI_Send(contactIdCountOfEachGrain.data(),numOfGrainsInRank,MPI_INT,0,rank,CART_COMM);
        MPI_Send(contactInfoCountOfEachGrain.data(),numOfGrainsInRank,MPI_INT,0,rank,CART_COMM);
      }
      MPI_Wait(&request, MPI_STATUS_IGNORE);
      /*send contact id of each grain in this rank*/
      int totalGrainContactIds = allContactPairOfRank.size();
      MPI_Isend(&totalGrainContactIds,1,MPI_INT,0,rank,CART_COMM,&request);
      if(totalGrainContactIds != 0){
        MPI_Send(allContactPairOfRank.data(),totalGrainContactIds,MPI_INT,0,rank,CART_COMM);
      }
      MPI_Wait(&request, MPI_STATUS_IGNORE);
    }else{
      // for rank 0, initialize coordinateNumberList
      for(int i = 0; i < num_grains; ++i){
        if(belongToThisRank[i] == 1){
          coordinateNumberList[i] = grainContactInfo[i].size();
          for(int otherId : grainContactId[i]){
            grainContactId[otherId].push_back(i);
          }
        }
      }
      for(int _rank = 1; _rank < USED; ++_rank){
        /*#################### FOR NORMAL GRAINS #####################*/
        std::vector<int> recvGrainIdsOfRank;
        std::vector<int> recvContactIdCountOfEachGrain;
        std::vector<int> recvContactInfoCountOfEachGrain;
        std::vector<int> RecvAllContactPairOfRank;
        MPI_Status status;
        int numOfGrainsInRank = 0;
        MPI_Recv(&numOfGrainsInRank,1,MPI_INT,_rank,_rank,CART_COMM,&status);
        if(numOfGrainsInRank != 0){
          recvGrainIdsOfRank.resize(numOfGrainsInRank);
          recvContactIdCountOfEachGrain.resize(numOfGrainsInRank);
          recvContactInfoCountOfEachGrain.resize(numOfGrainsInRank);
          MPI_Recv(recvGrainIdsOfRank.data(),numOfGrainsInRank,MPI_INT,_rank,_rank,CART_COMM,&status);
          MPI_Recv(recvContactIdCountOfEachGrain.data(),numOfGrainsInRank,MPI_INT,_rank,_rank,CART_COMM,&status);
          MPI_Recv(recvContactInfoCountOfEachGrain.data(),numOfGrainsInRank,MPI_INT,_rank,_rank,CART_COMM,&status);
        }
        int totalGrainContactIds = 0;
        MPI_Recv(&totalGrainContactIds,1,MPI_INT,_rank,_rank,CART_COMM,&status);
        if(totalGrainContactIds != 0){
          RecvAllContactPairOfRank.resize(totalGrainContactIds);
          MPI_Recv(RecvAllContactPairOfRank.data(),totalGrainContactIds,MPI_INT,_rank,_rank,CART_COMM,&status);
        }
        /*put in grainContactId*/
        int current = 0;
        for(int i = 0; i < numOfGrainsInRank; ++i){
          int gid = recvGrainIdsOfRank[i];
          originRankList[gid] = _rank;
          coordinateNumberList[gid] = recvContactInfoCountOfEachGrain[i];
          for(int j = 0; j < recvContactIdCountOfEachGrain[i]; ++j){
            int otherId = RecvAllContactPairOfRank[current + j];
            grainContactId[gid].push_back(otherId);
            grainContactId[otherId].push_back(gid);
          }
          current += recvContactIdCountOfEachGrain[i];
        }
      }
    }
  }
  /*##########################################################*/
  /*############## The CutHill-McKee Algorithm ###############*/
  /*##########################################################*/
  if(rank == 0){
    std::vector<int> degreeList(num_grains,0);
    std::vector<bool> checkList(num_grains,false);
    std::priority_queue<DegreeQueue, std::vector<DegreeQueue>, CompareDegree> iterateQueue;
    /*#################### PREPARE DEGREE List #####################*/
    for(int i = 0; i < num_grains; ++i){
      degreeList[i] = grainContactId[i].size();
      iterateQueue.emplace(i, grainContactId[i].size());
    }
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
  }
  auto end_time3 = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff3 = end_time3 - start_time3;
  double seconds3 = diff3.count();
  t3 += seconds3;

  /*##########################################################*/
  /*################## BROADCAST ISLANDS #####################*/
  /*##########################################################*/
  MPI_Barrier(CART_COMM);
  auto start_time4 = std::chrono::steady_clock::now();
  std::vector<int> islandRank;  // which island should be handled by which rank
  std::vector<int> eachIslandSize; // each island size
  std::vector<vector<int>> vectorizedSet;
  std::vector<int> allVectorizedSet;
  {
    //broadcast islands info, need islandsCount, count in each island, island rank
    if(rank == 0){
      islandRank.resize(islandsCount);
      vectorizedSet.resize(islandsCount);
      //determine which island should be handled by which rank
      for(int i = 0; i < islandsCount; ++i){
        std::vector<int> rankCount(USED, 0);
        for(auto gid : islands[i]){
          rankCount[originRankList[gid]] += coordinateNumberList[gid];
        }
        //find the maximum count
        int maxCount  = 0;
        int maxRankId = 0;
        for(int _rank = 0; _rank < USED; ++_rank){
          if(rankCount[_rank] > maxCount){
            maxRankId = _rank;
            maxCount  = rankCount[_rank];
          }
        }
        //decide the rankId
        islandRank[i] = maxRankId;
        eachIslandSize.push_back(islands[i].size());
        vectorizedSet[i].assign(islands[i].begin(), islands[i].end());
        allVectorizedSet.insert(allVectorizedSet.end(), islands[i].begin(), islands[i].end());
      }
    }
    //begin broadcast
    //islandsCount
    MPI_Bcast(&islandsCount,1,MPI_INT,0,CART_COMM);
    if(rank != 0){
      islandRank.resize(islandsCount);
      eachIslandSize.resize(islandsCount);
      vectorizedSet.resize(islandsCount);
      allVectorizedSet.resize(num_grains);
    }
    MPI_Bcast(islandRank.data(),islandsCount,MPI_INT,0,CART_COMM);
    MPI_Bcast(eachIslandSize.data(),islandsCount,MPI_INT,0,CART_COMM);
    MPI_Bcast(allVectorizedSet.data(),num_grains,MPI_INT,0,CART_COMM);
    //assign vectorizedSet in rank != 0
    if(rank != 0){
      int current = 0;
      for(int i = 0; i < islandsCount; ++i){
        vectorizedSet[i].assign(allVectorizedSet.begin()+current, allVectorizedSet.begin()+current+eachIslandSize[i]);
        current += eachIslandSize[i];
      }
    }
  }
  auto end_time4 = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff4 = end_time4 - start_time4;
  double seconds4 = diff4.count();
  t4 += seconds4;

  /*##########################################################*/
  /*################## ALLTOALL CONTACTS #####################*/
  /*##########################################################*/
  MPI_Barrier(CART_COMM);
  auto start_time5 = std::chrono::steady_clock::now();
  {
    std::vector<int> sendCount(USED,0);
    std::vector<int> sendDisp(USED,0);
    std::vector<int> recvCount(USED,0);
    std::vector<int> recvDisp(USED,0);
    std::vector<std::vector<ContactInfo>> sendGrainContactInfo(USED, std::vector<ContactInfo>());
    std::vector<ContactInfo> allSendGrainContactInfo;
    std::vector<ContactInfo> allRecvGrainContactInfo;
    // iterate sets, for each set, if current rank does not handle it
    // go inside the set, for each gid which is in current rank, send to buffer
    for(int i = 0; i < islandsCount; ++i){
      if(islandRank[i] != rank){
        for(auto gid : vectorizedSet[i]){
          if(belongToThisRank[gid] == 1){
            int rankId = islandRank[i];
            sendCount[rankId] += grainContactInfo[gid].size();
            sendGrainContactInfo[rankId].insert(sendGrainContactInfo[rankId].end(), grainContactInfo[gid].begin(), grainContactInfo[gid].end());
          }
        }
      }
    }
    //construct allSendGrainContactInfo
    for(int i = 0; i < USED; ++i){
      allSendGrainContactInfo.insert(allSendGrainContactInfo.end(),sendGrainContactInfo[i].begin(),sendGrainContactInfo[i].end());
    }
    MPI_Alltoall(sendCount.data(),1,MPI_INT,recvCount.data(),1,MPI_INT,CART_COMM);
    int totalRecv = 0;
    for(int i = 0; i < USED; ++i){
      totalRecv += recvCount[i];
      if(i == 0){
        sendDisp[i] = 0;
        recvDisp[i] = 0;
      }else{
        sendDisp[i] = sendDisp[i-1] + sendCount[i-1];
        recvDisp[i] = recvDisp[i-1] + recvCount[i-1];
      }
    }
    allRecvGrainContactInfo.resize(totalRecv);
    MPI_Alltoallv(allSendGrainContactInfo.data(),sendCount.data(),sendDisp.data(),CONTACT,
                  allRecvGrainContactInfo.data(),recvCount.data(),recvDisp.data(),CONTACT,CART_COMM);

    for(int i = 0; i < totalRecv; ++i){
      int master = allRecvGrainContactInfo[i]._master;
      int slave = allRecvGrainContactInfo[i]._slave;
      Vector3d normal = allRecvGrainContactInfo[i]._normal;
      if(belongToThisRank[master] != 1){
        grainContactInfo[master].push_back(allRecvGrainContactInfo[i]);
      }
    }
  }
  auto end_time5 = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff5 = end_time5 - start_time5;
  double seconds5 = diff5.count();
  t5 += seconds5;

  /*##########################################################*/
  /*################# COLLISION RESOLUTION ###################*/
  /*##########################################################*/
  //expensive, all deep copy...
  MPI_Barrier(CART_COMM);
  auto start_time6 = std::chrono::steady_clock::now();
  std::vector<vector<ContactInfo>> contactLists;
  for(int i = 0; i < islandsCount; ++i){
    if(islandRank[i] == rank){
      std::vector<ContactInfo> cList;
      for(auto & gid : vectorizedSet[i]){
        if(belongToThisRank[gid] != 2){
          cList.insert(cList.end(), grainContactInfo[gid].begin(), grainContactInfo[gid].end());
        }
      }
      contactLists.push_back(cList);
    }
  }

  int largestList = 0;
  int totalList = 0;
  for(int i = 0; i < contactLists.size(); ++i){
    largestList = largestList > contactLists[i].size() ? largestList : contactLists[i].size();
    totalList += contactLists[i].size();
  }
  MPI_Allreduce(MPI_IN_PLACE, &largestList, 1, MPI_INT, MPI_MAX, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &totalList, 1, MPI_INT, MPI_SUM, CART_COMM);
  if(rank == 0) cout<<"Largest contact list is: "<<largestList<<" #islands: "<<islandsCount<<" total contact: "<<totalList<<endl;


  /*##########################################################*/
  /*################# RIGID WALL ITERATION ###################*/
  /*##########################################################*/
  int tempN = 5;
  for(int tempStep = 0; tempStep < tempN; tempStep++){
    for(auto & contactList : contactLists){
      //collision_resolution(contactList,step);
      
      for(int j = 0; j < contactList.size() / 250; ++j){
        auto first = contactList.begin() + j * 250;
        auto last = contactList.begin() + (j+1) * 250;
        vector<ContactInfo> newVec(first, last);
        if(newVec.size() > 0) collision_resolution(newVec, step);
      }
      auto first = contactList.begin() + contactList.size() / 250 * 250;
      auto last = contactList.end();
      vector<ContactInfo> newVec(first, last);
      if(newVec.size() > 0) collision_resolution(newVec, step);

    }

    /* broadcast velocity to grain */
    MPI_Barrier(CART_COMM);
    auto start_time7 = std::chrono::steady_clock::now();
    {
      std::vector<int> recvCount(USED);
      std::vector<int> recvDisp(USED,0);
      std::vector<vel_data> velSendBuffer;
      std::vector<vel_data> velRecvBuffer(num_grains);
      //construct allGatherCount
      for(int i = 0; i < islandsCount; ++i){
        int rankId = islandRank[i];
        recvCount[rankId] += vectorizedSet[i].size();
      }
      //construct allGatherDisp
      for(int i = 1; i < USED; ++i){
        recvDisp[i] = recvDisp[i-1] + recvCount[i-1];
      }
      //construct allGatherSendBuffer
      for(int i = 0; i < islandsCount; ++i){
        if(islandRank[i] == rank){
          for(auto gid : vectorizedSet[i]){
            velSendBuffer.emplace_back(gid, rank, grainsWorld[gid].getVelocity(), grainsWorld[gid].getOmegaGlobal());
          }
        }
      }
      MPI_Allgatherv(velSendBuffer.data(),recvCount[rank],VELOCITY,
        velRecvBuffer.data(),recvCount.data(),recvDisp.data(),VELOCITY,CART_COMM);
      //update velocities
      for(auto vel : velRecvBuffer){
        if(vel._rank != rank){
          int gid = vel._id;
          grainsWorld[gid].changeVelocity(vel._velocity);
          grainsWorld[gid].changeOmegaGlobal(vel._omegaGlobal);
        }
      }
    }
    auto end_time7 = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff7 = end_time7 - start_time7;
    double seconds7 = diff7.count();
    t7 += seconds7;

    /* take temp time step */
    for(int i = 0; i < num_grains; ++i){
      if(belongToThisRank[i] == 1){
        grainsWorld[i].takeTempTimestep(dt);
      }
    }
    /* update grains positions, correct velocities */
    MPI_Barrier(CART_COMM);
    start_time7 = std::chrono::steady_clock::now();
    {
      set<int> grainList;
      vector<vel_data> boundaryVelSendBuffer;
      vector<vel_data> boundaryVelRecvBuffer;
      int sendCount = 0;
      int recvTotal = 0;
      vector<int> recvCount(USED);
      vector<int> recvDisp(USED,0);
      for(int i = 0; i < num_grains; ++i){
        if(belongToThisRank[i] == 1){
          if(wallPlanes[LEFT].bTempCircleCheck(grainsWorld[i]) ){
            wallPlanes[LEFT].findGrainWallContact(grainsWorld[i], dt);
            grainList.insert(i);
          }
          if(wallPlanes[BACK].bTempCircleCheck(grainsWorld[i]) ){
            wallPlanes[BACK].findGrainWallContact(grainsWorld[i], dt);
            grainList.insert(i);
          }
          if(wallPlanes[DOWN].bTempCircleCheck(grainsWorld[i]) ){
            wallPlanes[DOWN].findGrainWallContact(grainsWorld[i], dt);
            grainList.insert(i);
          }
          if(wallPlanes[RIGHT].bTempCircleCheck(grainsWorld[i]) ){
            wallPlanes[RIGHT].findGrainWallContact(grainsWorld[i], dt);
            grainList.insert(i);
          }
          if(wallPlanes[FRONT].bTempCircleCheck(grainsWorld[i]) ){
            wallPlanes[FRONT].findGrainWallContact(grainsWorld[i], dt);
            grainList.insert(i);
          }
          if(wallPlanes[UP].bTempCircleCheck(grainsWorld[i]) ){
            wallPlanes[UP].findGrainWallContact(grainsWorld[i], dt);
            grainList.insert(i);
          }
        }
      }
      for(auto & i : grainList){
        boundaryVelSendBuffer.emplace_back(i, rank, grainsWorld[i].getVelocity(), grainsWorld[i].getOmegaGlobal());
      }
      /* broadcast boundary velocity again */
      sendCount = boundaryVelSendBuffer.size();
      MPI_Allgather(&sendCount, 1, MPI_INT, recvCount.data(), 1, MPI_INT, CART_COMM);
      for(int i = 1; i < USED; ++i){
        recvDisp[i] = recvDisp[i-1] + recvCount[i-1];
      }
      for(int i = 0; i < USED; ++i){
        recvTotal += recvCount[i];
      }
      boundaryVelRecvBuffer.resize(recvTotal);
      MPI_Allgatherv(boundaryVelSendBuffer.data(), sendCount, VELOCITY,
          boundaryVelRecvBuffer.data(), recvCount.data(), recvDisp.data(), VELOCITY, CART_COMM);
      for(auto & vel : boundaryVelRecvBuffer){
        int gid = vel._id;
        if(belongToThisRank[gid] != 1){
          grainsWorld[gid].changeVelocity(vel._velocity);
          grainsWorld[gid].changeOmegaGlobal(vel._omegaGlobal);
        }
      }
    }
    end_time7 = std::chrono::steady_clock::now();
    diff7 = end_time7 - start_time7;
    seconds7 = diff7.count();
    t7 += seconds7;
    /* reset relative velocity */
    for(int i = 0; i < contactLists.size(); ++i){
      for(auto & item : contactLists[i]){
        item.computeRelVelocityThreshold();
      }
    }
    /* END LOOP */
  }
  auto end_time6 = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff6 = end_time6 - start_time6;
  double seconds6 = diff6.count();
  t6 += seconds6;
  MPI_Barrier(CART_COMM);
  auto start_time8 = std::chrono::steady_clock::now();
  /*##########################################################*/
  /*#################### POSITION UPDATE #####################*/
  /*##########################################################*/
  //takeTimestep
  for(int i = 0; i < num_grains; ++i){
    //all grain should take time step, maybe not?
    find_pos(grainsWorld[i].getPosition(), _idx, _idy, _idz);
    _old_bin_id = to_bin_id(_idx, _idy, _idz);
    grainsWorld[i].takeTimestep(dt);
    if(belongToThisRank[i]==1){
      find_pos(grainsWorld[i].getPosition(), _idx, _idy, _idz);
      _proc_id = to_proc_id(_idx, _idy, _idz);
      _bin_id = to_bin_id(_idx, _idy, _idz);
      if(_proc_id == rank){
        if(_bin_id != _old_bin_id){
          //cout<<"grain "<<i<<" move to different bin at rank: "<<rank<<" step: "<<step<<endl;
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
      }
      else{
        //enum DIRECTIONS {LEFT, RIGHT, FRONT, BACK, DOWN, UP};
        data_t migrated_grain = data_t(i,grainsWorld[i].getQuat(),
          grainsWorld[i].getPosition(),grainsWorld[i].getVelocity(),grainsWorld[i].getOmega());
        if(_proc_id%dimX < rank%dimX){
          grain_to_be_send[LEFT].push_back(migrated_grain);
          //cout<<"particle migration id: "<<i<<" step: "<<step<<" send to left. from rank: "<<rank<<endl;
        }else if(_proc_id%dimX > rank%dimX){
          grain_to_be_send[RIGHT].push_back(migrated_grain);
          //cout<<"particle migration id: "<<i<<" step: "<<step<<" send to right. from rank: "<<rank<<endl;
        }else if(_proc_id%(dimX*dimY)/dimX < rank%(dimX*dimY)/dimX){
          grain_to_be_send[FRONT].push_back(migrated_grain);
          //cout<<"particle migration id: "<<i<<" step: "<<step<<" send to front. from rank: "<<rank<<endl;
        }else if(_proc_id%(dimX*dimY)/dimX > rank%(dimX*dimY)/dimX){
          grain_to_be_send[BACK].push_back(migrated_grain);
          //cout<<"particle migration id: "<<i<<" step: "<<step<<" send to back. from rank: "<<rank<<endl;
        }else if(_proc_id/(dimX*dimY) < rank/(dimX*dimY)){
          grain_to_be_send[DOWN].push_back(migrated_grain);
          //cout<<"particle migration id: "<<i<<" step: "<<step<<" send to down. from rank: "<<rank<<endl;
        }else if(_proc_id/(dimX*dimY) > rank/(dimX*dimY)){
          grain_to_be_send[UP].push_back(migrated_grain);
          //cout<<"particle migration id: "<<i<<" step: "<<step<<" send to up. from rank: "<<rank<<endl;
        }

        //DELETE
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
    else if(belongToThisRank[i] == 2){
      // order here is important
      grainIdList[i] = -1;
      binIdList[i] = -1;
      belongToThisRank[i] = 0;
    }
  }

  updateAllPoints();
  auto end_time8 = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff8 = end_time8 - start_time8;
  double seconds8 = diff8.count();
  t8 += seconds8;
  //exchange grains
  MPI_Barrier(CART_COMM);
  auto start_time9 = std::chrono::steady_clock::now();
  if(USED > 1){
    grainsExchange(rank, indir);
  }
  auto end_time9 = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff9 = end_time9 - start_time9;
  double seconds9 = diff9.count();
  t9 += seconds9;

  MPI_Allreduce(MPI_IN_PLACE, &t1, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &t2, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &t3, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &t4, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &t5, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &t6, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &t7, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &t8, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &t9, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
  T1 += t1;
  T2 += t2;
  T3 += t3;
  T4 += t4;
  T5 += t5;
  T6 += t6;
  T7 += t7;
  T8 += t8;
  T9 += t9;

  //compute kinematic energy
  double kinematicEnergy = 0.0;
  int maxID = 0;
  double maxVelocity = 0.0;
  double maxOmg = 0.0;
  for(int i = 0; i < num_grains; ++i){
    if(belongToThisRank[i] == 1){
      kinematicEnergy += grainsWorld[i].computeKineticEnergy();
      maxID = maxVelocity < grainsWorld[i].getVelocity().norm() ? i : maxID;
      maxVelocity = maxVelocity < grainsWorld[i].getVelocity().norm() ? grainsWorld[i].getVelocity().norm() : maxVelocity;
      maxOmg = maxOmg < grainsWorld[i].getOmega().norm() ? grainsWorld[i].getOmega().norm() : maxOmg;
    }
  }

  // MPI calls  sendbuf       recvbuff          count     type       op       comm
  MPI_Allreduce(MPI_IN_PLACE, &kinematicEnergy, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &maxVelocity, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
  MPI_Allreduce(MPI_IN_PLACE, &maxOmg, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);

  if(step % 5 == 0){
    if(rank == 0){
      cout<<" T1 = " << T1 << " seconds. T2 = " << T2 << " seconds. T3 = "<<T3<<endl;
      cout<<" T4 = " << T4 << " seconds. T5 = " << T5 << " seconds. T6 = "<<T6<<endl;
      cout<<" T7 = " << T7 << " seconds. T8 = " << T8 << " seconds. T9 = "<<T9<<endl;
      cout<<"THIS IS ITERATION: "<<step<<" kinematic energy is: "<<kinematicEnergy<<" max velocity is: "<<maxVelocity<<" max Omega is: "<<maxOmg<<endl;
    }
  }
  }
}
#endif /*MPI_SOLVER_RIGID_HPP_*/
