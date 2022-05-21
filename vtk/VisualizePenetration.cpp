/*
 * VisualizeMembraneInitState.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: reid
 *  Modified on Jun 4, 2020
 *      Author: Peng TAN (UC Berkeley)
 */

#include "definitions.h"
#include "readInputMethods.h"
#include "buildMethods.h"
#include <vtkPlaneSource.h>
#include  <vtkSphereSource.h>
#include <random>
#include "omp.h"

int main(int argc, char *argv[])
{
	size_t ng = 394;
	string testName = "impulse_400";
	string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/InitState";
	vector<Vector3d> positions = readPositionFile(result_path+"/positions_"+testName+".dat", ng);
	vector<Vector4d> rotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", ng);

	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(1.0, 1.0, 1.0); // white
	renderWindow->AddRenderer(renderer);

	{
		Vector3d pos(0.,0.,45.);
		size_t gid = 182;
		Vector4d rot = rotations[gid];
		stringstream fname;
		fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<gid<<".dat";
		Polyhedron poly = readPolyFile2(fname.str());
		vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, pos, rot/rot.norm());
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(trianglePolyData);
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
		actor->GetProperty()->SetOpacity(1.);
		renderer->AddActor(actor);
	}

	{
		Vector3d pos(0.,0.,10);
		size_t   gid = 50;
		Vector4d rot = rotations[gid];
		stringstream fname;
		fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<gid<<".dat";
		Polyhedron poly = readPolyFile2(fname.str());
		vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, pos, rot/rot.norm());
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(trianglePolyData);
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(1., 0., 0.);
		actor->GetProperty()->SetOpacity(1.);
		renderer->AddActor(actor);
	}

/*
	{
		Vector3d pos(0.,0.,25.);
		Vector4d rot = rotations[gid];
		stringstream fname;
		fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<gid<<".dat";
		Polyhedron poly = readPolyFile2(fname.str());
		vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, pos, rot/rot.norm());
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(trianglePolyData);
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(1.0, 1.0, 0.0);
		actor->GetProperty()->SetOpacity(0.2);
		renderer->AddActor(actor);
	}

	{
		Vector3d pos(0.,10.,36.);
		Vector4d rot = rotations[gid-1];
		stringstream fname;
		fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<gid<<".dat";
		Polyhedron poly = readPolyFile2(fname.str());
		vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, pos, rot/rot.norm());
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(trianglePolyData);
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(1.0, 0., 0.);
		actor->GetProperty()->SetOpacity(1.);
		renderer->AddActor(actor);
	}
*/

	{
		vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
		planeSource->SetPoint1(200., 0., 0.);
		planeSource->SetPoint2(0., 200., 0.);
		planeSource->SetCenter(0., 0., 0.);
		planeSource->Update();
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(planeSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.15, 0., 0.);
		actor->GetProperty()->SetOpacity(0.25);
		renderer->AddActor(actor);
	}

	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
	renderer->SetActiveCamera(camera);
	camera->SetViewUp(0,0,1);
	camera->SetPosition(-7.5,0.,1.5);
	camera->SetFocalPoint(0,0,2.);
	camera->SetClippingRange(3,100);
	camera->Zoom(5.);
	renderWindow->Render();
	renderer->ResetCamera();

	//Screenshot
	vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
	windowToImageFilter->SetInput(renderWindow);
	windowToImageFilter->SetScale(3,3);
	windowToImageFilter->ReadFrontBufferOff();
	windowToImageFilter->Update();
	stringstream fname; fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/demo.png";

	//writer for saving later
	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(fname.str().c_str());
	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
	writer->Write();

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Start();

	cout << "program terminated" << endl;
	return EXIT_SUCCESS;
}
