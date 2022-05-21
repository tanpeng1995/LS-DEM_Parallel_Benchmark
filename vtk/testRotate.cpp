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
#include  <vtkSphereSource.h>
#include  <vtkCylinderSource.h>
#include  <vtkTransform.h>
#include  <vtkPlaneSource.h>
#include  <random>

namespace{
  vtkSmartPointer<vtkPolyData> CreateArrow(double& length, array<double, 3>& p1, array<double, 3>& p2);
}

int main(int argc, char *argv[])
{
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(1.0, 1.0, 1.0); // white
  renderWindow->AddRenderer(renderer);

  Vector3d position(10.,10.,10.);
  Vector3d force(100.,100.,141.4);

  //cylinder
  vtkNew<vtkCylinderSource> cylinder;
  cylinder->SetHeight(20.0);                           // height
  cylinder->SetRadius(1.0);     // radius
  cylinder->SetResolution(100);                        // number of slides
  cylinder->SetCenter(position(0), position(1), position(2));  // set center
  vtkNew<vtkPolyDataMapper> cylinderMapper;
  cylinderMapper->SetInputConnection(cylinder->GetOutputPort());
  vtkNew<vtkActor> clyinderActor;
  clyinderActor->SetMapper(cylinderMapper);
  clyinderActor->GetProperty()->SetColor(0.0,0.9,0.9);
  clyinderActor->GetProperty()->SetOpacity(1.);
  clyinderActor->SetOrigin(position(0), position(1), position(2));
  double alpha = asin(force(2)/force.norm());
  double beta  = atan(force(0)/force(1));
  cout<<"beta = "<<beta*180./3.1415926<<" alpha = "<<alpha*180.0/3.1415926<<endl;
  clyinderActor->RotateZ(-beta*180./3.1415926);
  if((force(0) < 0) && (force(1) < 0) || (force(0) > 0) && (force(1) < 0)) alpha = -alpha;
  clyinderActor->RotateX(alpha*180./3.1415926);

  //plane
  vtkNew<vtkPlaneSource> plane;
  plane->SetOrigin(0.,0.,0.);
  plane->SetPoint1(30.,0.,0.);
  plane->SetPoint2(0.,20.,0.);
  plane->Update();
  vtkNew<vtkPolyDataMapper> planeMapper;
  planeMapper->SetInputConnection(plane->GetOutputPort());
  vtkNew<vtkActor> planeActor;
  planeActor->SetMapper(planeMapper);
  planeActor->GetProperty()->SetColor(1.0,0.,0.);

  renderer->AddActor(clyinderActor);
  renderer->AddActor(planeActor);


  // Camera
  vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
  renderer->SetActiveCamera(camera);
  //Full View
  //camera->SetViewUp(0,1,0);
  /*
  camera->SetPosition(200,-100,275);
  camera->SetFocalPoint(25,25,250);
  camera->SetClippingRange(3,100);
  */
  // THIS FUCKING THING RIGHT HERE (mak sure that the object is between
  //the first value and the second value wrt the camera's position/focal point)
  camera->Zoom(.9);

  // Render and interact
  renderWindow->Render();
  renderer->ResetCamera();

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Start();

	cout << "program terminated" << endl;
	return EXIT_SUCCESS;
}
