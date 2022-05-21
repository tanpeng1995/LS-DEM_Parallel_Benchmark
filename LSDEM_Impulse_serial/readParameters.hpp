#ifndef READ_PARAMETERS_HPP_
#define READ_PARAMETERS_HPP_
#include "definitions.hpp"
#include "utilities.hpp"
#include "Grain3d.hpp"

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

std::unique_ptr<Grain3d> generateGrainFromFile(const Vector3d & position, const Vector4d & quat, const string & morphfile, double density, int grainidx) {
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
	// update grain object in the vector that was created at the beginning of this function
  //            mass   position   velocity   momentInertia   quat    omega   cmLset   refPointList  radius  lset mu  id
	return std::make_unique<Grain3d>(mass, density, position, velocity, momentOfInertia, quat, omega, cmLset, refPointList, radius, lset, mu, grainidx);
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

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  restart = atoi(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  tmax = atoi(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  tConsolidate = atoi(parameters_line.c_str());

  getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  dt = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  gDamping = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
  getline(parametersFile, parameters_line);
  threshold = atof(parameters_line.c_str());

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
  mu = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	wallFriction = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	wallStiffness = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	dHeight = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	ballRadius = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	pressure = atof(parameters_line.c_str());
	pressure = pressure * scalingFactor;

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	ballDensity = atof(parameters_line.c_str());
	ballDensity = ballDensity * scalingFactor * scalingFactor * scalingFactor;

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	kn = atof(parameters_line.c_str());

	getline(parametersFile, parameters_line,':');
	getline(parametersFile, parameters_line);
	kmn = atof(parameters_line.c_str());

	parametersFile.close();
}
#endif /* READ_PARAMETERS_HPP_ */
