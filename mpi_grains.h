/*
 * MPI_GRAINS_H_
 * Created on: July 13, 2021
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_GRAINS_H_
#define MPI_GRAINS_H_
#include "mpi_helper_function.h"
#include "definitions.h"
#include "utilities.h"

inline bool ghost_bin(int idx, int idy, int idz, int rank){
  int pz = rank/(dimX*dimY);
  int py = rank%(dimX*dimY)/dimX;
  int px = rank%(dimX*dimY)%dimX;

  return (idz == pz*(binsInDimZ-2)-1 && idy >= py*(binsInDimY-2)-1 && idy <= (py+1)*(binsInDimY-2)
          && idx >= px*(binsInDimX-2)-1 && idx <= (px+1)*(binsInDimX-2)) || /*bottom*/
         (idz == (pz+1)*(binsInDimZ-2) && idy >= py*(binsInDimY-2)-1 && idy <= (py+1)*(binsInDimY-2)
          && idx >= px*(binsInDimX-2)-1 && idx <= (px+1)*(binsInDimX-2)) || /*upper*/
         (idz >= pz*(binsInDimZ-2)-1 && idz <= (pz+1)*(binsInDimZ-2) && idy == py*(binsInDimY-2)-1
          && idx >= px*(binsInDimX-2)-1 && idx <= (px+1)*(binsInDimX-2)) || /*front*/
         (idz >= pz*(binsInDimZ-2)-1 && idz <= (pz+1)*(binsInDimZ-2) && idy == (py+1)*(binsInDimY-2)
          && idx >= px*(binsInDimX-2)-1 && idx <= (px+1)*(binsInDimX-2)) || /*back*/
         (idz >= pz*(binsInDimZ-2)-1 && idz <= (pz+1)*(binsInDimZ-2) && idy >= py*(binsInDimY-2)-1
          && idy <= (py+1)*(binsInDimY-2) && idx == px*(binsInDimX-2)-1) || /*left*/
         (idz >= pz*(binsInDimZ-2)-1 && idz <= (pz+1)*(binsInDimZ-2) && idy >= py*(binsInDimY-2)-1
          && idy <= (py+1)*(binsInDimY-2) && idx == px*(binsInDimX-2)-1) ; /*right*/
}

Grain3d regenerateGrainFromFile(const Vector3d & position, const Vector4d & quat, const Vector3d & velocity, const Vector3d & omega, const string & morphfile, int grainidx) {
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

	// create level set object
	Levelset3d lset = Levelset3d(lsetvec, xdim, ydim, zdim); //need copy constructor here
	// update grain object in the vector that was created at the beginning of this function
  //            mass   position   velocity   momentInertia   quat    omega   cmLset   refPointList  radius  lset  kn  ks  mu  id
	Grain3d grain(mass, position, velocity, momentOfInertia, quat, omega, cmLset, refPointList, radius, lset, kn, ks, mu, grainidx);
	return grain;
}

void filterUnusedGrain(){
  for(int i = 0; i < num_grains; ++i){
    if(grainsWorld[i].getRealGrain() && belongToThisRank[i] == 0){
      grainsWorld[i].clearGrain();
    }
  }
}
#endif /*MPI_GRAINS_H_*/
