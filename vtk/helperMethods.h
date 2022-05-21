/*
 * helperMethods.h
 *
 *  Created on: Nov 15, 2016
 *      Author: reid
 */

#ifndef HELPERMETHODS_H_
#define HELPERMETHODS_H_

#include "definitions.h"
#define PI 3.1415926;


double quaternionsDiff(Vector4d q1, Vector4d q2){
	q2 << q2(0), -q2(1), -q2(2), -q2(3);
	double a = q1(0)*q2(0) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3);
	double b = q1(0)*q2(1) + q1(1)*q2(0) + q1(2)*q2(3) - q1(3)*q2(2);
	double c = q1(0)*q2(2) - q1(1)*q2(3) + q1(2)*q2(0) + q1(3)*q2(1);
	double d = q1(0)*q2(3) + q1(1)*q2(2) - q1(2)*q2(1) + q1(3)*q2(0);
	double alpha = acos(a)*2;
	double beta1 = acos(b/sin(alpha/2));
	double beta2 = acos(c/sin(alpha/2));
	double beta3 = acos(d/sin(alpha/2));
	return alpha*180.0/PI;
}

// quaternion to rotation matrix
Vector4d qfromR(const Matrix3d & R) {
	Vector4d _quat;
	double tr = R.trace();
	_quat(3) = sqrt(tr+1.)/2.;
	_quat(0) = (R(0,2)-R(2,0))/(4*_quat(3));
	_quat(1) = (R(1,2)-R(2,1))/(4*_quat(3));
	_quat(2) = (R(0,1)-R(1,0))/(4*_quat(3));
	return _quat;
}

// rotation matrix to quaternion
Matrix3d Rfromq(Vector4d quat) { // pass by value
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

// finds a color from a rgb colormap linearly given some value and its max/min values
Vector3d findColor(const double & val, const double & minVal, const double & maxVal, const vector<Vector3d> & cmap) {
	Vector3d col;
	if (val <= minVal) {
		col = cmap[0];
	}
	else if (val >= maxVal) {
		col = cmap[cmap.size()-1];
	}
	else {
		double pctile = (val-minVal)/(maxVal-minVal);
		double idx = floor((double)cmap.size()*pctile);
		col = cmap[(size_t)idx];
	}
	return col;
}

// output camera information on keypress
void KeypressCallbackFunction (vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) ) {
	vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
	vtkRenderWindow *wren = iren->GetRenderWindow();
	vtkRendererCollection *rens = wren->GetRenderers();
	vtkRenderer *ren = rens->GetFirstRenderer();
	vtkCamera *cam = ren->GetActiveCamera();
	double pos[3]; double fpt[3]; double vup[3];
	cam->GetPosition(pos);
	cam->GetFocalPoint(fpt);
	cam->GetViewUp(vup);
	cout << "Camera position: " << pos[0] << ", " << pos[1] << ", " << pos[2] <<
			  " Focal point: " << fpt[0] << ", " << fpt[1] << ", " << fpt[2] <<
			  " View up: " << vup[0] << ", " << vup[1] << ", " << vup[2] << endl;
}


#endif /* HELPERMETHODS_H_ */
