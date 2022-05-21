/*
 * VisualizeFabric.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: reid
 *  Modified on Jun 4, 2020
 *      Author: Peng TAN (UC Berkeley)
 */

#include "definitions.h"
#include "readInputMethods.h"
#include "buildMethods.h"
#include  <vtkSphereSource.h>
#include <random>

int main(int argc, char *argv[])
{
  size_t gid = atoi(argv[1]);
	size_t ng = 569;
	string testName = "fabric_400";
	string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/InitState";
	vector<Vector3d> positions = readPositionFile(result_path+"/positions_"+testName+".dat", ng);
	vector<Vector4d> rotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", ng);

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  // render window interacter
  // create renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(1.0, 1.0, 1.0); // white
  renderWindow->AddRenderer(renderer);

  std::mt19937 gen(123123);
	std::uniform_real_distribution<double> rand_real(0., 1.0);

  stringstream fname;
  fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<gid+1<<".dat";
  Polyhedron poly = readPolyFile(fname.str());
  vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, positions[gid], rotations[gid]);
  // Create mapper and actor
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(trianglePolyData);

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  //actor->GetProperty()->SetColor(rand_real(gen), rand_real(gen), rand_real(gen));
  actor->GetProperty()->SetColor(0.8, 0.8, 0.8);
  actor->GetProperty()->SetOpacity(1);
  renderer->AddActor(actor);

  // Camera
	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
	renderer->SetActiveCamera(camera);
	//Full View
	camera->SetViewUp(0,0,-1);
	camera->SetPosition(-50,50,-50);
	camera->SetFocalPoint(20,0,0);
	camera->SetClippingRange(3,100);
	// THIS FUCKING THING RIGHT HERE (mak sure that the object is between
	//the first value and the second value wrt the camera's position/focal point)
	camera->Zoom(.9);

	// Render and interact
	renderWindow->Render();
	renderer->ResetCamera();

	//Screenshot
	vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
	windowToImageFilter->SetInput(renderWindow);
	windowToImageFilter->SetScale(3,3);
	windowToImageFilter->ReadFrontBufferOff();
	windowToImageFilter->Update();
	stringstream fname2;
  fname2 << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/Fabric_"<<gid<<".png";

	//writer for saving later
  vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(fname2.str().c_str());
	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
	writer->Write();

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Start();

	cout << "program terminated" << endl;
	return EXIT_SUCCESS;
}
