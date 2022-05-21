/*
 * MainTest.cpp
 *
 * Created on: May 6, 2020
 *      Author: Peng TAN (Berkeley)
 */
#include "definitions.h"
#include "utilities.h"
#include "mpi_gather_for_save.h"
#include "mpi_solver.h"
#include "mpi_init_simulation.h"

void readParametersTrueUnit(string parametersdir){
  ifstream parametersFile(parametersdir);
	string parameters_filepath;
	string parameters_line;
	if(parametersFile.fail()){
		cout<<"FAIL TO OPEN FILE: Parameters_LSDEM.txt"<<endl;
		return;
	}
	while (parameters_line.compare("## INPUT PARAMETERS BEGIN HERE ##")){
    getline(parametersFile, parameters_line);
  }

  getline(parametersFile, parameters_line,'"');
  getline(parametersFile, testName,'"');

	getline(parametersFile, parameters_line,'"');
  getline(parametersFile, outdirName,'"');

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  restart = atoi(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  massScaling = atoi(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  tmax = atoi(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  tConsolidate = atoi(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  timeFac = atof(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  outputFrequency = atoi(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  scalingFactor = atof(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  num_grains = atoi(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  grainDensity = atof(parameters_line.c_str());
	grainDensity = grainDensity*scalingFactor*scalingFactor*scalingFactor;

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  kn = atof(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  ks = atof(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  mu = atof(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  gDamping = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  resolution = atoi(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  worldWidth = atof(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  worldLength = atof(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  worldHeight = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  dynamic_binning = atoi(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  binning_frequency = atoi(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  initHeight = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  initRadius = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  ballRadius = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	pressure = atof(parameters_line.c_str());
	pressure = pressure*scalingFactor;

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	wallDensity = atof(parameters_line.c_str());
	wallDensity = wallDensity*scalingFactor*scalingFactor*scalingFactor;

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	wallToGrainStiffness = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	kmn = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	kms = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	wallFriction = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	capHeight = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	capRadius = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	capDensity = atof(parameters_line.c_str());
	capDensity = capDensity*scalingFactor*scalingFactor*scalingFactor;

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	dHeight = atof(parameters_line.c_str());

	parametersFile.close();
}

int main(int argc, char* argv[]) {
  string parametersdir = "/global/home/users/tanpeng/LSDEM_Parallel_ball/Parameters_LSDEM.txt";
  readParametersTrueUnit(parametersdir);
  string indir = "/global/scratch/users/tanpeng/LSDEM_Modelling/Input/"+testName+"/";
  string outdir = "/global/scratch/users/tanpeng/LSDEM_Modelling/Output/"+outdirName+"/";

  FILE * posfile  = fopen((outdir+"positions_"+testName+".dat").c_str(), "w");
  FILE * rotfile  = fopen((outdir+"radius_"+testName+".dat").c_str(), "w");
  FILE * wallfile = fopen((outdir+"walls_"+testName+".dat").c_str(), "w");
  FILE * wallfile2 = fopen((outdir+"external_walls_"+testName+".dat").c_str(), "w");
  FILE * stressFile = fopen((outdir+"stress_"+testName+".dat").c_str(), "w");
  FILE * volumeCNFile = fopen((outdir+"volume_coordination_"+testName+".dat").c_str(), "w");

  init_simulation(indir, wallfile, wallfile2);

  auto start_time = std::chrono::steady_clock::now();
  const double dt = timeFac*sqrt(minMass/kn);

  cout<<"################## LS-DEM-BALL #####################"<<endl;
  cout<<"Begin Simulation and Good Luck, Author: Peng TAN (Berkeley)"<<endl;
  cout<<"massScaling: "<<massScaling<<" minMass: "<<minMass<<" maxMass: "<<maxMass<<" dt: "<<dt<<endl;
  cout<<"This code is serial, it models grains as balls, used for comparison"<<endl;
  cout<<"For strain control simulation, dHeight = "<<dHeight<<endl;
  cout<<"Using different friction coefficient, mu: "<<mu<<endl;
  cout<<"Output Directory: "<<outdir<<endl;

  for (int step = 0; step < tmax; ++step) {
    //STEP 1: Simulate One Step
    simulate_one_step(dt, step, indir);
    //STEP 2: Bookkeeping
    if(step % outputFrequency == 0){
      gather_for_save(step, outdir, indir, posfile, rotfile, wallfile, wallfile2, stressFile, volumeCNFile);
    }
  }

  auto end_time = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end_time - start_time;
  double seconds = diff.count();

  cout<<"Simulation Time = " << int(floor(seconds/3600.)) << " hrs, "<< int(floor(fmod(seconds/60., 60.))) << " mins, "<< -floor(fmod(seconds/60., 60.))*60. + fmod(seconds, 3600.) << " seconds\n";

  fclose(posfile);
  fclose(rotfile);
  fclose(wallfile);
  fclose(wallfile2);
  fclose(volumeCNFile);

  return 0;
}
