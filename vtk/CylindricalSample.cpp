/*
 *  try to obtain an accurate cylindrical shape from a sample
 *  Created on: Nov 17, 2021
 *      Author: Peng TAN (Berkeley)
 */
#include "definitions.h"
#include <fstream>

// quaternion to rotation matrix
Vector4d RotationMatToQuat(const Matrix3d & R) {
	Vector4d _quat;
	double tr = R.trace();
	_quat(3) = sqrt(tr+1.)/2.;
	_quat(0) = (R(0,2)-R(2,0))/(4*_quat(3));
	_quat(1) = (R(1,2)-R(2,1))/(4*_quat(3));
	_quat(2) = (R(0,1)-R(1,0))/(4*_quat(3));
	return _quat;
}

// rotation matrix to quaternion
Matrix3d QuatToRotationMat(Vector4d quat) { // pass by value
	Matrix3d _rotMatrix;
	_rotMatrix << 1., 0., 0., 0., 1., 0., 0., 0., 1.;
	if(isinf(quat(0)) || isinf((-quat(0))) ) return _rotMatrix;
	quat = quat/quat.norm();
	_rotMatrix(0,0) = -quat(0)*quat(0) + quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
	_rotMatrix(0,1) = -2*(quat(0)*quat(1) - quat(2)*quat(3));
	_rotMatrix(0,2) =  2*(quat(1)*quat(2) + quat(0)*quat(3));
	_rotMatrix(1,0) = -2*(quat(0)*quat(1) + quat(2)*quat(3));
	_rotMatrix(1,1) =  quat(0)*quat(0) - quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
	_rotMatrix(1,2) = -2*(quat(0)*quat(2) - quat(1)*quat(3));
	_rotMatrix(2,0) =  2*(quat(1)*quat(2) - quat(0)*quat(3));
	_rotMatrix(2,1) = -2*(quat(0)*quat(2) + quat(1)*quat(3));
	_rotMatrix(2,2) = -quat(0)*quat(0) - quat(1)*quat(1) + quat(2)*quat(2) + quat(3)*quat(3);
	return _rotMatrix;
}

vector<Vector3d> getPositions(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string  line;
	string 	partial;
	istringstream iss;
	vector<Vector3d> positions;
	positions.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
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

vector<Vector4d> getQuaternions(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<Vector4d> quaternions;
	quaternions.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
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

bool qualifiedGrain(const Vector3d & position, const Vector4d & quat,
	const string & morphfile, int gid, const double & lowHeight,
	const double & highHeight, const double & radius, const double & minMass, const double & maxMass){
	ifstream file(morphfile.c_str());
	if(file.fail()){
		cout << "Grain "<< gid + 1 << " does not exist." << endl;
		return false;
	}
	string line;
	string partial;
	istringstream iss;
	//mass
	getline(file, line);
	double mass = atof(line.c_str());
	if(mass < minMass || mass > maxMass) {
		return false;
	}
	// moment of inertia
	getline(file, line);
	//cmLset
	getline(file, line);
	// npoints
	getline(file, line);
	int npoints = atoi(line.c_str());
	// points
	Matrix3d R = QuatToRotationMat(quat);
	getline(file, line);
	Vector3d point;
	iss.str(line);
	for(int i = 0; i < npoints; ++i){
		getline(iss, partial, ' ');
		point(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(2) = atof(partial.c_str());
		point = R * point + position;
		if(point(2) < lowHeight || point(2) > highHeight || point(0)*point(0) + point(1)*point(1) > radius*radius){
			return false;
		}
	}
	iss.clear();
	return true;
}


bool qualifiedGrain2(const Vector3d & position, const Vector4d & quat,
	const string & morphfile, int gid, const double & lowHeight,
	const double & highHeight, const double & width, const double & length, const double & minMass, const double & maxMass){
	ifstream file(morphfile.c_str());
	if(file.fail()){
		cout << "Grain "<< gid + 1 << " does not exist." << endl;
		return false;
	}
	string line;
	string partial;
	istringstream iss;
	//mass
	getline(file, line);
	double mass = atof(line.c_str());
	if(mass < minMass || mass > maxMass) {
		return false;
	}
	// moment of inertia
	getline(file, line);
	//cmLset
	getline(file, line);
	// npoints
	getline(file, line);
	int npoints = atoi(line.c_str());
	// points
	Matrix3d R = QuatToRotationMat(quat);
	getline(file, line);
	Vector3d point;
	iss.str(line);
	for(int i = 0; i < npoints; ++i){
		getline(iss, partial, ' ');
		point(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(2) = atof(partial.c_str());
		point = R * point + position;
		if(point(2) < lowHeight || point(2) > highHeight || abs(point(0)) > width || abs(point(1)) > length ){
			return false;
		}
	}
	iss.clear();
	return true;
}

int main2(int argc, char* argv[]){
	size_t ng = atoi(argv[1]);
  string testName = "7-6B";
  string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/InitState";

	vector<Vector3d> positions = getPositions(result_path+"/positions_"+testName+".dat", ng);
	vector<Vector4d> rotations = getQuaternions(result_path+"/rotations_"+testName+".dat", ng);

	const double lowHeight  = 100.0;
	const double highHeight = 1100.0;
	const double width      = 400.0;
	const double length     = 1000.0;
	const double minMass    = 400.0;
	const double maxMass    = 20000.0;

	vector<int> qualifiedId;
	for(int i = 0; i < ng; ++i){
		stringstream fname;
		if(i % 50 == 0) cout<<"have proceed "<<i<<" grains."<<endl;
		fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Morphologies/morph_"<<i+1<<".dat";
		if( qualifiedGrain2(positions[i], rotations[i], fname.str(), i, lowHeight, highHeight, width, length, minMass, maxMass) ){
			qualifiedId.push_back(i);
		}
	}

	//save
	std::ofstream outFile("/home/hasitha/Desktop/data/fabric/7-6B_1000/qualified_id.txt");
  // the important part
  for (const auto &id : qualifiedId) outFile << id << " ";
	return 0;
}


int main(int argc, char* argv[]){
	size_t ng = atoi(argv[1]);
  string testName = "full_belgium";
  string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/InitState";

	vector<Vector3d> positions = getPositions(result_path+"/positions_"+testName+".dat", ng);
	vector<Vector4d> rotations = getQuaternions(result_path+"/rotations_"+testName+".dat", ng);

	const double lowHeight  = 0.0;
	const double highHeight = 750.0;
	const double radius     = 375.0;
	const double minMass    = 3000.0;
	const double maxMass    = 150000.0;

	vector<int> qualifiedId;
	for(int i = 0; i < ng; ++i){
		stringstream fname;
		fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Morphologies/morph_"<<i+1<<".dat";
		if( qualifiedGrain(positions[i], rotations[i], fname.str(), i, lowHeight, highHeight, radius, minMass, maxMass) ){
			qualifiedId.push_back(i);
		}else{
			cout<<"grain "<<i<<" is not qualified."<<endl;
		}
	}

	//save
	std::ofstream outFile("/home/hasitha/Desktop/data/fabric/cylinder_belgium/qualified_id.txt");
  // the important part
  for (const auto &id : qualifiedId) outFile << id << " ";
	return 0;
}
