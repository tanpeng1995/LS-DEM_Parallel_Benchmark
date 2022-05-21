/*
 * readInputMethods.h
 *
 *  Created on: Nov 15, 2016
 *      Author: reid
 */

#include "definitions.h"

#ifndef READINPUTMETHODS_H_
#define READINPUTMETHODS_H_

// generally these objects will be in the reference config with centroid at 0
struct Polyhedron {
	Polyhedron() {}
	Polyhedron(vector<Vector3i> faces, vector<Vector3d> verts, vector<Vector3d> normals):
		_faces(faces), _verts(verts), _normals(normals) {}
	vector<Vector3i> _faces;
	vector<Vector3d> _verts;
	vector<Vector3d> _normals;
};

// reads file.  File format:
// nvertices
// v1x v1y v1z v2x v2y v2z ... vnx vny vnz [vertex coordinates]
// nfaces
// f11 f12 f13 f21 f22 f23 ... fn1 fn2 fn3 [faces]
// n1x n1y n1z n2x n2y n2z ... nnx nny nnz [vertex normals]
Polyhedron readPolyFile(string filename) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	// vertices
	getline(file, line);
	size_t nverts = atoi(line.c_str());
	vector<Vector3d> verts(nverts);
	for (size_t n = 0; n < nverts; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		verts[n](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		verts[n](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		verts[n](2) = atof(partial.c_str());
	}
	// faces
	getline(file, line);
	size_t nfaces = atoi(line.c_str());
	vector<Vector3i> faces(nfaces);
	for (size_t n = 0; n < nfaces; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		faces[n](0) = atoi(partial.c_str());
		getline(iss, partial, ' ');
		faces[n](1) = atoi(partial.c_str());
		getline(iss, partial, ' ');
		faces[n](2) = atoi(partial.c_str());
	}
	// vertex normals
	vector<Vector3d> normals(nverts);
	for (size_t n = 0; n < nverts; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		normals[n](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normals[n](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normals[n](2) = atof(partial.c_str());
		normals[n] = normals[n]/normals[n].norm();
	}
	Polyhedron poly(faces, verts, normals);
	return poly;
}

Polyhedron readTopographyFile(string filename) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	// vertices
	getline(file, line);
	size_t nverts = atoi(line.c_str());
	vector<Vector3d> verts(nverts);
	for (size_t n = 0; n < nverts; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		verts[n](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		verts[n](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		verts[n](2) = atof(partial.c_str());
	}
	// faces
	getline(file, line);
	size_t nfaces = atoi(line.c_str());
	vector<Vector3i> faces(nfaces);
	for (size_t n = 0; n < nfaces; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		faces[n](0) = atoi(partial.c_str());
		getline(iss, partial, ' ');
		faces[n](1) = atoi(partial.c_str());
		getline(iss, partial, ' ');
		faces[n](2) = atoi(partial.c_str());
	}
	// vertex normals
	vector<Vector3d> normals(nverts);
	for (size_t n = 0; n < nverts; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		normals[n](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normals[n](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normals[n](2) = atof(partial.c_str());
		normals[n] = normals[n]/normals[n].norm();
	}
	Polyhedron poly(faces, verts, normals);
	cout<<"Topography data: #faces = "<<faces.size()<<" #verts = "<<verts.size()<<" #normals = "<<normals.size()<<endl;
	return poly;
}


Polyhedron readPolyFile2(string filename) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	// mass
	getline(file, line);
	// vertices
	getline(file, line);
	size_t nverts = atoi(line.c_str());
	vector<Vector3d> verts(nverts);
	for (size_t n = 0; n < nverts; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		verts[n](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		verts[n](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		verts[n](2) = atof(partial.c_str());
	}
	// faces
	getline(file, line);
	size_t nfaces = atoi(line.c_str());
	vector<Vector3i> faces(nfaces);
	for (size_t n = 0; n < nfaces; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		faces[n](0) = atoi(partial.c_str());
		getline(iss, partial, ' ');
		faces[n](1) = atoi(partial.c_str());
		getline(iss, partial, ' ');
		faces[n](2) = atoi(partial.c_str());
	}
	// vertex normals
	vector<Vector3d> normals(nverts);
	for (size_t n = 0; n < nverts; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		normals[n](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normals[n](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normals[n](2) = atof(partial.c_str());
		normals[n] = normals[n]/normals[n].norm();
	}
	Polyhedron poly(faces,verts, normals);
	return poly;
}

Polyhedron readPolyFile3(string filename, double & mass) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	// mass
	getline(file, line);
	mass = atof(line.c_str());
	// vertices
	getline(file, line);
	size_t nverts = atoi(line.c_str());
	vector<Vector3d> verts(nverts);
	for (size_t n = 0; n < nverts; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		verts[n](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		verts[n](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		verts[n](2) = atof(partial.c_str());
	}
	// faces
	getline(file, line);
	size_t nfaces = atoi(line.c_str());
	vector<Vector3i> faces(nfaces);
	for (size_t n = 0; n < nfaces; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		faces[n](0) = atoi(partial.c_str());
		getline(iss, partial, ' ');
		faces[n](1) = atoi(partial.c_str());
		getline(iss, partial, ' ');
		faces[n](2) = atoi(partial.c_str());
	}
	// vertex normals
	vector<Vector3d> normals(nverts);
	for (size_t n = 0; n < nverts; n++) {
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		normals[n](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normals[n](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normals[n](2) = atof(partial.c_str());
		normals[n] = normals[n]/normals[n].norm();
	}
	Polyhedron poly(faces,verts, normals);
	return poly;
}

vector<Vector3d> readPositionFile(string filename, size_t ngrains) {
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

vector<Vector3d> readPositionFile(string filename) {
	ifstream file(filename.c_str());
	string  line;
	string 	partial;
	istringstream iss;
	getline(file, line);
	size_t ngrains = atoi(line.c_str());
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
} // end

vector<Vector3d> readVectorFile(string filename) {
	ifstream file(filename.c_str());
	string  line;
	string 	partial;
	istringstream iss;
	getline(file, line);
	size_t N = atoi(line.c_str());
	vector<Vector3d> ret(N);
	for (size_t i = 0; i < N; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		ret[i](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		ret[i](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		ret[i](2) = atof(partial.c_str());
		iss.clear();
	}
	return ret;
}


vector<std::tuple<Vector3d, Vector3d>> readVector3dTupleFile(string filename){
	ifstream file(filename.c_str());
	string  line;
	string 	partial;
	istringstream iss;
	getline(file, line);
	size_t N = atoi(line.c_str());
	vector<std::tuple<Vector3d, Vector3d>> ret;
	for (size_t i = 0; i < N; i++) {
		Vector3d pos1(0.,0.,0.), pos2(0.,0.,0.);
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		pos1(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		pos1(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		pos1(2) = atof(partial.c_str());

		getline(iss, partial, ' ');
		pos2(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		pos2(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		pos2(2) = atof(partial.c_str());
		iss.clear();

		ret.push_back(std::make_tuple(pos1, pos2));
	}
	return ret;
}

vector<Vector4d> readQuaternionFile(string filename, size_t ngrains) {
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

vector<int> readIntegerFile(string filename, size_t num) {
	ifstream file(filename.c_str());
	string   line;
	vector<int> integers(num);
	for (size_t n = 0; n < num; n++) {
		getline(file, line);
		integers[n]= atoi(line.c_str());
	}
	return integers;
}

vector<double> readDoubleFile(string filename, size_t num) {
	ifstream file(filename.c_str());
	string   line;
	vector<double> doubles(num);
	for (size_t n = 0; n < num; n++) {
		getline(file, line);
		doubles[n]= atof(line.c_str());
	}
	return doubles;
}

vector<std::set<int>> readSetFile(string filename, size_t rank){
	vector<std::set<int>> retSet(rank);
	ifstream file(filename.c_str());
	string line;
	string partial;
	istringstream iss;
	getline(file, line);
	size_t num = atoi(line.c_str());
	cout<<"total #ball contacts: "<<num<<endl;
	int ballId, rankId;
	for (size_t i = 0; i < num; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		ballId = atoi(partial.c_str());
		getline(iss, partial, ' ');
		rankId = atoi(partial.c_str());
		iss.clear();
		retSet[rankId].insert(ballId);
	}
	return retSet;
}

struct VGposition{
	// VG = variable number of grains
	size_t ngrains;
	vector<Vector3d> positions;
	vector<size_t> idlist;

};

vector<VGposition> readVGpositionfile(string filename, size_t nstations) {
	ifstream file(filename.c_str());
	string line;
	string partial;
	istringstream iss;
	vector<VGposition> VGpos(nstations);
	Vector3d pos;
	for (size_t t=0; t<nstations;t++){
		getline(file,line);
		VGpos[t].ngrains = atoi(line.c_str());
		for (size_t n=0;n<VGpos[t].ngrains;n++){
			getline(file,line);
			iss.str(line);
			getline(iss, partial, ' ');
			pos(0) = atof(partial.c_str());
			getline(iss, partial, ' ');
			pos(1) = atof(partial.c_str());
			getline(iss, partial, ' ');
			pos(2) = atof(partial.c_str());
			iss.clear();
			VGpos[t].positions.push_back(pos);
		}
	}
	return VGpos;
}

struct VGrotation{
	// VG = variable number of grains
	size_t ngrains;
	vector<Vector4d> rotations;
	vector<size_t> idlist;

};

vector<VGrotation> readVGrotationfile(string filename, size_t nstations) {
	ifstream file(filename.c_str());
	string line;
	string partial;
	istringstream iss;
	vector<VGrotation> VGrot(nstations);
	Vector4d rot;
	for (size_t t=0; t<nstations;t++){
		getline(file,line);
		VGrot[t].ngrains = atoi(line.c_str());
		for (size_t n=0;n<VGrot[t].ngrains;n++){
			getline(file,line);
			iss.str(line);
			getline(iss, partial, ' ');
			rot(0) = atof(partial.c_str());
			getline(iss, partial, ' ');
			rot(1) = atof(partial.c_str());
			getline(iss, partial, ' ');
			rot(2) = atof(partial.c_str());
			getline(iss, partial, ' ');
			rot(3) = atof(partial.c_str());
			iss.clear();
			VGrot[t].rotations.push_back(rot);
		}
	}
	return VGrot;
}

void readIdListFile(string filename, vector<VGposition> & positions, size_t nstations){
	ifstream file(filename.c_str());
	string line;
	for (size_t t=0;t<nstations;t++){
		getline(file,line);
		for (size_t n=0;n<positions[t].ngrains;n++){
			getline(file,line);
			positions[t].idlist.push_back(atoi(line.c_str()));
		}
	}
}

struct CData {
	// member variables
	vector<size_t> _Ncontacts;
	vector<Vector2i> _cpairs;
	vector<Vector3d> _forces;
	vector<Vector3d> _normals;
	vector<Vector3d> _clocs;

	void operator+=(const CData & c) {
		_Ncontacts.insert( _Ncontacts.end(), c._Ncontacts.begin() , c._Ncontacts.end());
		_cpairs.insert( _cpairs.end(),	c._cpairs.begin(),	c._cpairs.end());
		_forces.insert( _forces.end(),	c._forces.begin(),	c._forces.end());
		_normals.insert(_normals.end(),	c._normals.begin(),	c._normals.end());
		_clocs.insert(_clocs.end(),		c._clocs.begin(),		c._clocs.end());
	}
};

CData readCdatafile(string filename, size_t nstations) {
	ifstream file(filename.c_str());
	string line;
	string partial;
	istringstream iss;
	CData cdata;
	Vector2i cpair;
	Vector3d force;
	Vector3d normal;
	Vector3d cloc;

	for (size_t t = 0;t<nstations;t++) {
		getline(file,line);
		cdata._Ncontacts.push_back( atoi(line.c_str()));
		for (size_t n = 0;n<cdata._Ncontacts[t];n++){
			//cpairs
			getline(file,line);
			iss.str(line);
			getline(iss, partial, ' ');
			cpair(0) = atoi(partial.c_str());
			getline(iss, partial, ' ');
			cpair(1) = atoi(partial.c_str());
			cdata._cpairs.push_back(cpair);
			iss.clear();
			//forces
			getline(file,line);
			iss.str(line);
			getline(iss, partial, ' ');
			force(0) = atoi(partial.c_str());
			getline(iss, partial, ' ');
			force(1) = atoi(partial.c_str());
			getline(iss, partial, ' ');
			force(2) = atoi(partial.c_str());
			cdata._forces.push_back(force);
			iss.clear();
			//normals
			getline(file,line);
			iss.str(line);
			getline(iss, partial, ' ');
			normal(0) = atoi(partial.c_str());
			getline(iss, partial, ' ');
			normal(1) = atoi(partial.c_str());
			getline(iss, partial, ' ');
			normal(2) = atoi(partial.c_str());
			cdata._normals.push_back(normal);
			iss.clear();
			//clocs
			getline(file,line);
			iss.str(line);
			getline(iss, partial, ' ');
			cloc(0) = atoi(partial.c_str());
			getline(iss, partial, ' ');
			cloc(1) = atoi(partial.c_str());
			getline(iss, partial, ' ');
			cloc(2) = atoi(partial.c_str());
			cdata._clocs.push_back(cloc);
			iss.clear();
		}
	}

	return cdata;

}

vector<Vector3d> readSphereMembrane(string filename, size_t & nb, double & ball_radius, size_t ns){
	ifstream file(filename.c_str());
	string line;
	string partial;
	istringstream iss;

	getline(file, line);
	nb = atoi(line.c_str());

	getline(file, line);
	ball_radius = atof(line.c_str());

	vector<Vector3d> ballPositions(nb*ns);
	for(int i = 0; i < nb*ns; ++i){
		iss.clear();
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		ballPositions[i](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		ballPositions[i](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		ballPositions[i](2) = atof(partial.c_str());
	}

	return ballPositions;
}

#endif /* READINPUTMETHODS_H_ */
