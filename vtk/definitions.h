/*
 * definitions.h
 *
 *  Created on: Nov 15, 2016
 *      Author: reid
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

// Eigen
#include <Eigen/Dense>
using Eigen::Vector3d; // position, velocity, etc
using Eigen::Vector2i; // for some stuff with 2 integers as elements
using Eigen::Vector3i; // 3-integer vector
using Eigen::Vector4d; // quaternions
using Eigen::Matrix3d; // rotation matrices
typedef Eigen::Matrix<double,1,6> Vector6d; // voigt stress
typedef Eigen::Matrix<int,1,6> Vector6i; // 6-integer vector

#include <tuple>
#include <map>
#include <set>

// OpenMP
#include <omp.h>

// C libraries
#include <float.h>		// DBL_MAX, INT_MAX
#include <stdio.h>		// printf, fgets
#include <stdlib.h>		// atoi
#include <omp.h>			// OpenMP

// C++ libraries
#include <cmath>
using std::sqrt;
using std::sin;
using std::cos;
using std::isnan;
using std::acos;
using std::asin;

#include <cstdlib>
using std::size_t;

#include <string>		// strings
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

#include<tr1/array>
using std::tr1::array;	// standard array


// VTK
#include <vtkVersion.h> // version
#include <vtkSmartPointer.h> // Smart pointers used in almost everything
#include <vtkCellData.h>
#include <vtkDoubleArray.h> // used by normals
#include <vtkTriangle.h> // used by faces
#include <vtkCellArray.h>
#include <vtkPolyData.h> // vtk polygon objects
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkAlgorithm.h>
#include <vtkPolygon.h>
#include <vtkLineSource.h>
// Display
#include <vtkCallbackCommand.h> // callbacks
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleJoystickCamera.h>
// Screenshot
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkPostScriptWriter.h>


#endif /* DEFINITIONS_H_ */
