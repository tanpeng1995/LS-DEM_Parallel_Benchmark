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

  string parametersdir = "/global/home/users/tanpeng/LSDEM_Membrane/Parameters_LSDEM.txt";
  readParametersTrueUnit(parametersdir);
  string indir = "/global/scratch/users/tanpeng/LSDEM_Modelling/Input/"+testName+"/";
  string outdir = "/global/scratch/users/tanpeng/LSDEM_Modelling/Output/"+outdirName+"/";

  FILE * posfile  = fopen((outdir+"positions_"+testName+".dat").c_str(), "w");
  FILE * rotfile  = fopen((outdir+"rotations_"+testName+".dat").c_str(), "w");
  FILE * wallfile = fopen((outdir+"walls_"+testName+".dat").c_str(), "w");
  FILE * wallfile2 = fopen((outdir+"external_walls_"+testName+".dat").c_str(), "w");
  FILE * stressFile = fopen((outdir+"stress_"+testName+".dat").c_str(), "w");
  FILE * volumeCNFile = fopen((outdir+"volume_coordination_"+testName+".dat").c_str(), "w");

  int rank, num_procs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  init_simulation(rank, num_procs, indir, wallfile, wallfile2);

  auto start_time = std::chrono::steady_clock::now();
  double T_force = 0.0;
  double T_border = 0.0;
  double T_migrate = 0.0;
  double T_comm = 0.0;
  const double dt = timeFac*sqrt(minMass/kn);

  if (rank == 0){
    cout<<"################## LS-DEM #####################"<<endl;
    cout<<"Begin Simulation and Good Luck, Author: Peng TAN (Berkeley)"<<endl;
    cout<<"wallToGrainStiffness: "<<wallToGrainStiffness<<" global damping: "<<gDamping<<" resolution: "<<resolution<<" coordinateNumber: "<<coordinateNumber<<endl;
    cout<<"max radius: "<<maxR<<" max mass: "<<maxMass<<" min mass: "<<minMass<<" massScaling: "<<massScaling<<endl;
    cout<<"Used processor(s) is: "<<USED<<" Time-Step: "<<dt<<" tmax: "<<tmax<<endl;
    cout<<"Model membrane with an extra processor"<<endl;
    cout<<"For strain control simulation, dHeight = "<<dHeight<<" Confining pressure is: "<<pressure/scalingFactor<<endl;
    cout<<"Output Directory: "<<outdir<<endl;
    cout<<"InitState read from: "<<indir+InitState<<endl;
    cout<<"loadingStages: "<<endl;
    for(auto item : loadingStages){ cout<<item<<" "; }
    double inertiaNumber = dHeight/dt/initHeight * avgRadius / sqrt(pressure/grainDensity);
    cout<<"average grain radius is: "<<avgRadius<<" inertiaNumber is: "<<inertiaNumber<<endl;
    cout<<endl;
    cout<<"################## LS-DEM CYCLIC TEST ##################"<<endl<<endl;;
  }

  for (int step = 0; step < tmax; ++step) {
    //STEP 1: Simulate One Step
    simulate_one_step(rank, dt, step, indir, outdir, volumeCNFile, T_force, T_border, T_migrate, T_comm);
    //STEP 2: Bookkeeping
    if(step % outputFrequency == 0){
      gather_for_save(rank, step, outdir, indir, posfile, rotfile, wallfile, wallfile2, stressFile);
    }
    //STEP 3: Dynamically Divide Current Computational Domain
    if(dynamic_binning && step > 0 && step % binning_frequency == 0 && num_procs > 1){
      if(rank == 0){ cout<<"dynamic_binning at step: "<<step<<endl; }
      reassign_particles(rank,indir);
    }
  }

  auto end_time = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end_time - start_time;
  double seconds = diff.count();
  if (rank == 0){
    cout<<"num_procs is: "<<num_procs<<" Simulation Time = " << int(floor(seconds/3600.)) << " hrs, "<< int(floor(fmod(seconds/60., 60.))) << " mins, "<< -floor(fmod(seconds/60., 60.))*60. + fmod(seconds, 3600.) << " seconds\n";
    cout<<"num_procs is: "<<num_procs<<" Force Interaction Time = " << int(floor(T_force/3600.)) << " hrs, "<< int(floor(fmod(T_force/60., 60.))) << " mins, "<< -floor(fmod(T_force/60., 60.))*60. + fmod(T_force, 3600.) << " seconds\n";
    cout<<"num_procs is: "<<num_procs<<" Border Communication Time = " << int(floor(T_border/3600.)) << " hrs, "<< int(floor(fmod(T_border/60., 60.))) << " mins, "<< -floor(fmod(T_border/60., 60.))*60. + fmod(T_border, 3600.) << " seconds\n";
    cout<<"num_procs is: "<<num_procs<<" Particle Migration Time = " << int(floor(T_migrate/3600.)) << " hrs, "<< int(floor(fmod(T_migrate/60., 60.))) << " mins, "<< -floor(fmod(T_migrate/60., 60.))*60. + fmod(T_migrate, 3600.) << " seconds\n";
  }

  fclose(posfile);
  fclose(rotfile);
  fclose(wallfile);
  fclose(wallfile2);
  fclose(volumeCNFile);

  if(rank < USED){
    freeMPIBuff_static();
    freeMPIBuff_dynamic();
  }

  MPI_Finalize();
  return 0;
}
