/*
 * MPI_BALLWALLS_CONTACT_H_
 * Created on: June 20, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_BALLWALLS_CONTACT_H_
#define MPI_BALLWALLS_CONTACT_H_

extern int openmpThreads;
extern int dynamicSchedule;

void collectBoundaryGrains(vector<boundary_data_send> & boundaryDataSend){
  boundaryDataSend.clear(); boundaryDataSend.shrink_to_fit();
  #pragma omp parallel for schedule(dynamic, dynamicSchedule) num_threads(openmpThreads)
  for(int i = 0; i < num_grains; ++i){
    if(belongToThisRank[i] == 1){
      #pragma omp critical
      {
        boundaryDataSend.emplace_back(i, grainsWorld[i].getQuat(), grainsWorld[i].getPosition());
      }
    }
  }
}

#endif /*MPI_BALLWALLS_CONTACT_H_*/
