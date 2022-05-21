/*
 * buildMethods.h
 *
 *  Created on: Nov 15, 2016
 *      Author: reid
 */

#include "definitions.h"
#include "helperMethods.h"
#include "readInputMethods.h"

#ifndef BUILDMETHODS_H_
#define BUILDMETHODS_H_

// creates a vtk polygon object from a Polygon object, first of many overloaded methods
vtkSmartPointer<vtkPolyData> constructVtkPoly(Polyhedron poly) {
	// vertices and vertex normals
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> normalsArray = vtkSmartPointer<vtkDoubleArray>::New();
	normalsArray->SetNumberOfComponents(3);
	normalsArray->SetNumberOfTuples(poly._verts.size());
	double vN[3];
	for (size_t v = 0; v < poly._verts.size(); v++) {
		points->InsertNextPoint (poly._verts[v](0), poly._verts[v](1), poly._verts[v](2));
		vN[0] = poly._normals[v](0); vN[1] = poly._normals[v](1); vN[2] = poly._normals[v](2);
		normalsArray->SetTuple(v, vN);
	}
	// faces
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	for (size_t f = 0; f < poly._faces.size(); f++) {
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		triangle->GetPointIds()->SetId ( 0, poly._faces[f](0) );
		triangle->GetPointIds()->SetId ( 1, poly._faces[f](1) );
		triangle->GetPointIds()->SetId ( 2, poly._faces[f](2) );
		triangles->InsertNextCell ( triangle );
	}
	// Create a polydata object
	vtkSmartPointer<vtkPolyData> polyhedron = vtkSmartPointer<vtkPolyData>::New();
	// Add the geometry and topology to the polydata
	polyhedron->SetPoints ( points );
	polyhedron->SetPolys ( triangles );
	polyhedron->GetPointData()->SetNormals(normalsArray);
	return polyhedron;
}

// creates a vtk polygon object from a Polygon object, a position, and a rotation
vtkSmartPointer<vtkPolyData> constructVtkPoly(Polyhedron poly, Vector3d pos, Vector4d rot) {
	// vertices and vertex normals
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> normalsArray = vtkSmartPointer<vtkDoubleArray>::New();
	normalsArray->SetNumberOfComponents(3);
	normalsArray->SetNumberOfTuples(poly._verts.size());
	double vN[3];
	Vector3d point;
	Vector3d normal;
	Matrix3d R = Rfromq(rot);
	for (size_t v = 0; v < poly._verts.size(); v++) {
		point = R*poly._verts[v] + pos;
		points->InsertNextPoint (point(0), point(1), point(2));
		normal = R*poly._normals[v];
		vN[0] = normal(0); vN[1] = normal(1); vN[2] = normal(2);
		normalsArray->SetTuple(v, vN);
	}
	// faces
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	for (size_t f = 0; f < poly._faces.size(); f++) {
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		triangle->GetPointIds()->SetId ( 0, poly._faces[f](0) );
		triangle->GetPointIds()->SetId ( 1, poly._faces[f](1) );
		triangle->GetPointIds()->SetId ( 2, poly._faces[f](2) );
		triangles->InsertNextCell ( triangle );
	}
	// Create a polydata object
	vtkSmartPointer<vtkPolyData> polyhedron = vtkSmartPointer<vtkPolyData>::New();
	// Add the geometry and topology to the polydata
	polyhedron->SetPoints ( points );
	polyhedron->SetPolys ( triangles );
	polyhedron->GetPointData()->SetNormals(normalsArray);
	return polyhedron;
}

// creates a vtk polygon object from a Polygon object, a position, and a rotation
vtkSmartPointer<vtkPolyData> constructVtkPoly(Polyhedron poly, const Vector3d & pos, const double & angle) {
	// vertices and vertex normals
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> normalsArray = vtkSmartPointer<vtkDoubleArray>::New();
	normalsArray->SetNumberOfComponents(3);
	normalsArray->SetNumberOfTuples(poly._verts.size());
	double vN[3];
	Vector3d point;
	Vector3d normal;
	Matrix3d R;
	R << cos(angle), 0., sin(angle),
	     0.,         1.,    0.,
			-sin(angle), 0., cos(angle);
	cout<<"rotation matrix of topography is: "<<endl;
	cout<<R<<endl;
	for (size_t v = 0; v < poly._verts.size(); v++) {
		point = R*poly._verts[v] + pos;
		points->InsertNextPoint (point(0), point(1), point(2));
		normal = R*poly._normals[v];
		vN[0] = normal(0); vN[1] = normal(1); vN[2] = normal(2);
		normalsArray->SetTuple(v, vN);
	}
	// faces
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	for (size_t f = 0; f < poly._faces.size(); f++) {
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		triangle->GetPointIds()->SetId ( 0, poly._faces[f](0) );
		triangle->GetPointIds()->SetId ( 1, poly._faces[f](1) );
		triangle->GetPointIds()->SetId ( 2, poly._faces[f](2) );
		triangles->InsertNextCell ( triangle );
	}
	// Create a polydata object
	vtkSmartPointer<vtkPolyData> polyhedron = vtkSmartPointer<vtkPolyData>::New();
	// Add the geometry and topology to the polydata
	polyhedron->SetPoints ( points );
	polyhedron->SetPolys ( triangles );
	polyhedron->GetPointData()->SetNormals(normalsArray);
	return polyhedron;
}


#endif /* BUILDMETHODS_H_ */
