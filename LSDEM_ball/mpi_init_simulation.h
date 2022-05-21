
/*
 * MPI_INIT_SIMULATION_H_
 *
 * Created on: June 20, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_INIT_SIMULATION_H_
#define MPI_INIT_SIMULATION_H_

vector<double> readDoubleFile(string filename, int num) {
	ifstream file(filename.c_str());
	string   line;
	vector<double> doubles(num);
	for (int n = 0; n < num; n++) {
		getline(file, line);
		doubles[n]= atof(line.c_str());
	}
	return doubles;
}
vector<Vector3d> readPositionFile(string filename, int ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	 partial;
	istringstream iss;
	vector<Vector3d> positions;
	positions.resize(ngrains);
	for (int grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		positions[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		positions[grainidx](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		positions[grainidx](2) = atof(partial.c_str());
		iss.clear();
	}
	return positions;
}

void init_simulation(string indir, FILE * wallfile, FILE * wallfile2){
  vector<Vector3d> inputPos;
  vector<double> inputRadius;
  inputPos = readPositionFile(indir + "InitState/positions_"+testName+".dat", num_grains);
  inputRadius = readDoubleFile(indir + "InitState/radius_"+testName+".dat", num_grains);
  grainsWorld = std::make_unique<Grain3d[]>(num_grains);
  cout<<"finish constructing linked-list..."<<endl;

  wallPlane = std::make_shared<WallPlane> (Vector3d(0.,0.,1.), Vector3d(0.,0.,0.), kn, wallFriction, initRadius, num_grains);
  wallCap = std::make_shared<WallCap>(Vector3d(0.,0.,-1.), Vector3d(0.,0.,initHeight), kn, wallFriction, num_grains+1, gDamping, capHeight, capRadius, capDensity);
	wallCylinder = std::make_shared<WallCylinder>(initHeight, initRadius, kn, ks, mu, num_grains+2);
  wallBalls = std::make_shared<WallBalls> (initHeight, initRadius, ballRadius, 1.0*pressure, wallDensity, wallToGrainStiffness, kmn, kms, gDamping);
  externalWall = std::make_shared<WallBalls>(initHeight, initRadius+2.0*ballRadius, 1.1*ballRadius, 0.2*pressure, wallDensity, wallToGrainStiffness, kmn, kms, gDamping);
  fprintf(wallfile, "%d\n", wallBalls->getNumWalls());
  fprintf(wallfile, "%.4f\n", wallBalls->getBallRadius());
  fprintf(wallfile2, "%d\n", externalWall->getNumWalls());
  fprintf(wallfile2, "%.4f\n", externalWall->getBallRadius());
  fflush(wallfile); fflush(wallfile2);
  cout<<"finish buiding wallPlane, wallCap, and wallBalls."<<endl;

  stressVoigt<<0., 0., 0., 0., 0., 0.;

  for(int i = 0; i < num_grains; ++i){
    if(i % 1000 == 0) cout<<"finish loading "<<i<<" grains"<<endl;
    grainsWorld[i] = Grain3d(inputPos[i], inputRadius[i], grainDensity, kn, ks, mu, i);
    minMass = min(minMass, grainsWorld[i].getMass());
    maxMass = max(maxMass, grainsWorld[i].getMass());
  }
  /* mass scaling */
  if(massScaling){
    minMass = maxMass / 20.;
    for(int i = 0; i < num_grains; ++i){
      if(grainsWorld[i].getMass() < minMass){
        grainsWorld[i].changeMass(minMass);
      }
    }
  }

}
#endif /*MPI_INIT_SIMULATION_H_*/
