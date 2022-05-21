/*
 * MainTest.cpp
 * Serial code for impulse-based algorithm
 * Mainly for debuggin
 * Created on: Feb 6, 2022
 *      Author: Peng TAN (Berkeley)
 */
#include "Levelset3d.hpp"
#include "Grain3d.hpp"
#include "WallBalls.hpp"
#include "WallPlane.hpp"
#include "utilities.hpp"
#include "collision_resolution.hpp"
#include "gather_for_save.hpp"
#include "init_simulation.hpp"
#include "solver.hpp"
#include "command_line_option.hpp"

int main(int argc, char* argv[]) {
  int openmpThreads = find_int_arg(argc, argv, "-n", 1);
  int dynamicSchedule = find_int_arg(argc, argv, "-d", 1);
  omp_set_num_threads(openmpThreads);

  string parametersdir = "/global/home/users/tanpeng/LSDEM_Impulse_serial/Parameters_LSDEM.txt";
  readParametersTrueUnit(parametersdir);
  string indir = "/global/scratch/users/tanpeng/LSDEM_Modelling/Input/"+testName+"/";
  string outdir = "/global/scratch/users/tanpeng/LSDEM_Modelling/Output/"+testName+"_2/";

  FILE * posfile  = fopen((outdir+"positions_"+testName+".dat").c_str(), "w");
  FILE * rotfile  = fopen((outdir+"rotations_"+testName+".dat").c_str(), "w");
  FILE * wallfile1 = fopen((outdir+"walls_1_"+testName+".dat").c_str(), "w");
  FILE * wallfile2 = fopen((outdir+"walls_2_"+testName+".dat").c_str(), "w");

  init_simulation(indir, wallfile1, wallfile2);

  auto start_time = std::chrono::steady_clock::now();
  cout<<"################## LS-DEM #####################"<<endl;
  cout<<"Begin Simulation and Good Luck, Author: Peng TAN (Berkeley)"<<endl;
  cout<<"Time-Step: "<<dt<<" tmax: "<<tmax<<endl;
  cout<<"Using impulse-based DEM, dHeight = "<<dHeight<<endl;
  cout<<"Output Directory: "<<outdir<<endl;

  for (int step = 0; step < tmax; ++step) {
    //STEP 1: Bookkeeping
    if(step % outputFrequency == 0){
      gather_for_save(step, outdir, posfile, rotfile, wallfile1, wallfile2);
    }
    //STEP 2: Simulate One Step
    simulate_one_step(dt, step, indir);
  }

  auto end_time = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end_time - start_time;
  double seconds = diff.count();

  cout<<" Simulation Time = " << int(floor(seconds/3600.)) << " hrs, "<< int(floor(fmod(seconds/60., 60.))) << " mins, "<< -floor(fmod(seconds/60., 60.))*60. + fmod(seconds, 3600.) << " seconds\n";

  fclose(posfile);
  fclose(rotfile);
  fclose(wallfile1);
  fclose(wallfile2);
  return 0;
}
