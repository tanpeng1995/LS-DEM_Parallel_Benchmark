#ifndef DEFINITIONS_HPP_
#define DEFINITIONS_HPP_
#include <chrono>
#include <random>
// C libraries
#include <float.h>		// DBL_MAX, INT_MAX
#include <stdio.h>		// printf, fgets
#include <stdlib.h>		// atoi
#include <omp.h>			// OpenMP
#include <mpi.h>			// MPI
#include <memory>
#include <memory>
#include <time.h>
#include <sys/time.h>
#include <queue>
// C++ libraries
#include <cmath>
using std::min;
using std::sqrt;
using std::sin;
using std::cos;
using std::isnan;
using std::abs;
using std::pow;
#include <cstdlib>
using std::size_t;
#include <cstring>		// strings
using std::string;
#include <iostream>		// terminal output
using std::cout;
using std::endl;
#include <fstream>		// file io
using std::stringstream;
using std::ifstream;
using std::getline;
using std::istringstream;
#include <algorithm>	// math stuff
using std::min;
using std::max;
#include <vector>		// standard vector
using std::vector;
#include <set>
using std::set;
#include <set>
using std::set;
#include<tr1/array>
using std::tr1::array;	// standard array
#include <Eigen/Dense>
#include <Eigen/SVD>
using Eigen::MatrixXd; // shellWalls
using Eigen::Vector3d; // position, velocity, etc
using Eigen::Vector2d;
using Eigen::Vector2i; // for some stuff with 2 integers as elements
using Eigen::Vector3i; // 3-integer vector
using Eigen::Vector4d; // quaternions
//typedef Eigen::Matrix<double, 4, 1, Eigen::DontAlign> Vector4d;
using Eigen::Matrix3d; // rotation matrices
using Eigen::Matrix2d;
//using Eigen::aligned_allocator;
//#include <Eigen/StdVector>
typedef Eigen::Matrix<double,1,6> Vector6d; // voigt stress
typedef Eigen::Matrix<int,1,6> Vector6i; // 6-integer vector
#endif /* DEFINITIONS_HPP_ */
