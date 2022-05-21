/*
 * MainTest.cpp
 *
 * Created on: May 6, 2020
 *      Author: Peng TAN (Berkeley)
 */
#include "mpi_helper_function.h"
#include "mpi_dynamic_binning.h"
#include "command_line_option.h"
#include "mpi_init_simulation.h"
#include "mpi_solver.h"
#include "mpi_gather_for_save.h"


int main(int argc, char* argv[]) {
  openmpThreads = find_int_arg(argc, argv, "-n", 1);
  dynamicSchedule = find_int_arg(argc, argv, "-d", 1);
  omp_set_num_threads(openmpThreads);

  string parametersdir = "/global/home/users/tanpeng/LSDEM_Test_duplicate/Parameters_LSDEM.txt";
  readParametersTrueUnit(parametersdir);
  string indir = "/global/scratch/users/tanpeng/LSDEM_Modelling/Input/"+testName+"/";
  string outdir = "/global/scratch/users/tanpeng/LSDEM_Modelling/Output/"+testName+"_6/";

  FILE * posfile  = fopen((outdir+"positions_"+testName+".dat").c_str(), "w");
  FILE * rotfile  = fopen((outdir+"rotations_"+testName+".dat").c_str(), "w");

  int rank, num_procs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  init_simulation(rank, num_procs, indir);

  auto start_time = std::chrono::steady_clock::now();
  double T_force  = 0.0;
  double T_border = 0.0;
  double T_migrate= 0.0;
  double T_update = 0.0;
  const double dt = timeFac*sqrt(minMass/kn);

  if (rank == 0){
    cout<<"################## LS-DEM #####################"<<endl;
    cout<<"Begin Simulation and Good Luck, Author: Peng TAN (Berkeley)"<<endl;
    cout<<"Used processor(s) is: "<<USED<<" Time-Step: "<<dt<<" tmax: "<<tmax<<endl;
    cout<<"Simuation type is settlement in cylinder"<<endl;
    cout<<"For strain control simulation, dHeight = "<<dHeight<<endl;
    cout<<"Output Directory: "<<outdir<<endl;
  }


  for (int step = 0; step < tmax; ++step) {
    //STEP 1: Bookkeeping
    if(step % outputFrequency == 0){
      gather_for_save(rank, posfile, rotfile, outdir, step);
    }
    //STEP 2: Simulate One Step
    simulate_one_step(rank, dt, step, indir, T_force, T_border, T_migrate, T_update);
  }

  auto end_time = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end_time - start_time;
  double seconds = diff.count();
  if (rank == 0){
    cout<<"num_procs is: "<<num_procs<<" Simulation Time = " << seconds << " seconds\n";
    cout<<"num_procs is: "<<num_procs<<" Force Interaction Time = " << T_force << " seconds\n";
    cout<<"num_procs is: "<<num_procs<<" Border Communication Time = " << T_border << " seconds\n";
    cout<<"num_procs is: "<<num_procs<<" Particle Migration Time = " << T_migrate << " seconds\n";
    cout<<"num_procs is: "<<num_procs<<" Info Update Time = " << T_update << " seconds\n";
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
