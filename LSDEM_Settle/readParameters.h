/*
 * READ_PARAMETERS_H_
 *
 *  Created on: May 30, 2014
 *      Author: Reid Kawamoto (Caltech)
 *  Modified on: May 6, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef READ_PARAMETERS_H_
#define READ_PARAMETERS_H_
#include "definitions.h"
#include "utilities.h"

vector<int> readIntegerFile(string filename){
	ifstream file(filename.c_str());
	string   line;
	vector<int> integers;
	while (getline(file, line))
	  integers.push_back(atoi(line.c_str()));
	return integers;
}

vector<int> readIntegerFile(string filename, int num) {
	ifstream file(filename.c_str());
	string   line;
	vector<int> integers(num);
	for (int n = 0; n < num; n++){
		getline(file, line);
		integers[n]= atoi(line.c_str());
	}
	return integers;
}

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

	vector<Vector3d> readVector3dFile(string filename, int ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<Vector3d> positions(ngrains);
	for (int grainidx = 0; grainidx < ngrains; ++grainidx) {
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

vector<Vector3d> readVector3dFile(string filename) {
	ifstream file(filename.c_str());
	string   line;
	string 	 partial;
	istringstream iss;
	vector<Vector3d> positions;
	int grainidx = 0;
	while (getline(file, line)) {
		iss.str(line);
		getline(iss, partial, ' ');
		positions.push_back( Vector3d(0,0,0) );
		positions[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		positions[grainidx](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		positions[grainidx](2) = atof(partial.c_str());
		iss.clear();
		grainidx++;
	}
	return positions;
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

vector<Vector4d> readQuaternionFile(string filename, int ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	 partial;
	istringstream iss;
	vector<Vector4d> quaternions;
	quaternions.resize(ngrains);
	for (int grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		quaternions[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		quaternions[grainidx](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		quaternions[grainidx](2) = atof(partial.c_str());
		getline(iss, partial, ' ');
		quaternions[grainidx](3) = atof(partial.c_str());
		iss.clear();
	}
	return quaternions;
}

// get properties.  order: density, kn, ks, mu
vector<double> readPropertiesFile(string filename) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<double> properties(4);
	properties.resize(4);
	for (int p = 0; p < 4; p++) {
		getline(file, line);
		properties[p]= atof(line.c_str());
	}
	return properties;
}

std::unique_ptr<Grain3d> generateGrainFromFile(const Vector3d & position, const Vector4d & quat, const string & morphfile, int grainidx) {
  ifstream file(morphfile.c_str());
	if(file.fail()) cout << "Grain " << grainidx+1 << " does not exist." << endl;
	string  line;
	string 	partial;
	istringstream iss;
	// temp stuff
	Vector3d momentOfInertia;
	Vector3d point;
	Vector3d cmLset;
	// mass
	getline(file, line);
	double mass = atof(line.c_str());
	// moment of inertia
	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	momentOfInertia(0) = atof(partial.c_str());
	getline(iss, partial, ' ');
	momentOfInertia(1) = atof(partial.c_str());
	getline(iss, partial, ' ');
	momentOfInertia(2) = atof(partial.c_str());
	iss.clear();
	// cmLset
	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	cmLset(0) = atof(partial.c_str());
	getline(iss, partial, ' ');
	cmLset(1) = atof(partial.c_str());
	getline(iss, partial, ' ');
	cmLset(2) = atof(partial.c_str());
	iss.clear();
	// npoints (INTEGER)
	getline(file, line);
	int npoints = atoi(line.c_str());
	// points
	getline(file, line);
	vector<Vector3d> refPointList(npoints);
	iss.str(line);
	for (int ptidx = 0; ptidx < npoints; ptidx++) {
		getline(iss, partial, ' ');
		point(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(2) = atof(partial.c_str());
		refPointList[ptidx] = point;
	}
	iss.clear();
	// bbox radius
	getline(file, line);
	double radius = atof(line.c_str());
	// level set dimensions (INTEGERS)
	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	int xdim = atoi(partial.c_str());
	getline(iss, partial, ' ');
	int ydim = atoi(partial.c_str());
	getline(iss, partial, ' ');
	int zdim = atoi(partial.c_str());
	iss.clear();
	// level set
	getline(file, line);
	vector<float> lsetvec(xdim*ydim*zdim);
	iss.str(line);
	for (int i = 0; i < xdim*ydim*zdim; i++) {
		getline(iss, partial, ' ');
		lsetvec[i] = atof(partial.c_str());
	}
	iss.clear();
	file.close();
	const Vector3d velocity(0.,0.,0.);
	const Vector3d omega(0.,0.,0.);
	// create level set object
	Levelset3d lset = Levelset3d(lsetvec, xdim, ydim, zdim); //need copy constructor here
  //             mass  position  velocity   momentInertia   quat  omega  cmLset  refPointList  radius  lset  kn  ks  mu  id
	return std::make_unique<Grain3d>(mass, position, velocity, momentOfInertia, quat, omega, cmLset, refPointList, radius, lset, kn, ks, mu, grainidx);
}

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
  massScaling = atoi(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  coordinateNumber = atoi(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  resolution = atoi(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  tmax = atoi(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  timeFac = atof(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  outputFrequency = atoi(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  dynamicBinning = atoi(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  binningFrequency = atoi(parameters_line.c_str());

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
  worldWidth = atof(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  worldLength = atof(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  worldHeight = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  gravity = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  initHeight = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  initRadius = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  endHeight = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  dHeight = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  wallFriction = atof(parameters_line.c_str());

	parametersFile.close();
}
#endif /* READ_PARAMETERS_H_ */
