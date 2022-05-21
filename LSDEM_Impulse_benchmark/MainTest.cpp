/*
 * MainTest.cpp
 *
 * Created on: May 6, 2020
 *      Author: Peng TAN (Berkeley)
 */
#include "mpi_helper_function.hpp"
#include "mpi_dynamic_binning.hpp"
#include "command_line_option.hpp"
#include "mpi_init_simulation.hpp"
#include "mpi_solver_rigid.hpp"
#include "mpi_gather_for_save.hpp"

//int openmpThreads = 4;
//int dynamicSchedule = 1;
int main(int argc, char* argv[]) {
  int openmpThreads = find_int_arg(argc, argv, "-n", 1);
  int dynamicSchedule = find_int_arg(argc, argv, "-d", 1);
  omp_set_num_threads(openmpThreads);

  string parametersdir = "/global/home/users/tanpeng/LSDEM_Impulse_benchmark/Parameters_LSDEM.txt";
  readParametersTrueUnit(parametersdir);
  string indir = "/global/scratch/users/tanpeng/LSDEM_Modelling/Input/"+testName+"/";
  string outdir = "/global/scratch/users/tanpeng/LSDEM_Modelling/Output/"+testName+"/";

  FILE * posfile  = fopen((outdir+"positions_"+testName+".dat").c_str(), "w");
  FILE * rotfile  = fopen((outdir+"rotations_"+testName+".dat").c_str(), "w");

  int rank, num_procs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  init_simulation(rank, num_procs, indir);

  auto start_time = std::chrono::steady_clock::now();
  double T1 = 0., T2 = 0., T3 = 0., T4 = 0., T5 = 0., T6 = 0., T7 =0., T8 = 0., T9 = 0.;

  if (rank == 0){
    cout<<"################## LS-DEM #####################"<<endl;
    cout<<"Begin Simulation and Good Luck, Author: Peng TAN (Berkeley)"<<endl;
    cout<<"Used processor(s) is: "<<USED<<" Time-Step: "<<dt<<" tmax: "<<tmax<<endl;
    cout<<"For mass scaling, minMass: "<<minMass<<" maxMass: "<<maxMass<<endl;
    cout<<"Using impulse-based DEM, parallel benchmark "<<endl;
    cout<<"Output Directory: "<<outdir<<endl;
  }

  for (int step = 0; step < tmax; ++step) {
    //STEP 1: Bookkeeping
    if(step % outputFrequency == 0){
      gather_for_save(rank, step, outdir, posfile, rotfile);
    }
    //STEP 2: Simulate One Step
    simulate_one_step(rank, dt, step, indir, T1, T2, T3, T4, T5, T6, T7, T8, T9);
    //STEP 3: Dynamically Divide Current Computational Domain
    if(dynamic_binning && step > 0 && step % binning_frequency == 0 && num_procs > 1){
      if(rank == 0){ cout<<"dynamic_binning at step: "<<step<<endl; }
      reassign_particles(rank,indir);
    }
    if(step == 500){
      T1 = 0.; T2 = 0.; T3 = 0.; T4 = 0.; T5 = 0.; T6 = 0.; T7 = 0.; T8 = 0.; T9 = 0.;	
    }
  }

  auto end_time = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end_time - start_time;
  double seconds = diff.count();

  if (rank == 0){
    cout<<"num_procs is: "<<num_procs<<" Simulation Time = " << int(floor(seconds/3600.)) << " hrs, "<< int(floor(fmod(seconds/60., 60.))) << " mins, "<< -floor(fmod(seconds/60., 60.))*60. + fmod(seconds, 3600.) << " seconds\n";
  }

  fclose(posfile);
  fclose(rotfile);

  if(rank < USED){
    freeMPIBuff_static();
    freeMPIBuff_dynamic();
  }

  MPI_Finalize();
  return 0;
}
