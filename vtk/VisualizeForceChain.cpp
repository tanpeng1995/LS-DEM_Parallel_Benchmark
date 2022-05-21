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
#include  <random>
#include <vtkRenderLargeImage.h>

int main(int argc, char *argv[])
{
  size_t ns = atoi(argv[1]);
  double mmf = atof(argv[2]);
  string testName = "7-6B_cylinder";
  string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/results";
  vector<Vector3d> forces    = readVectorFile(result_path+"/forceChain_"+testName+"_"+std::to_string(ns)+".dat");
	vector<std::tuple<Vector3d, Vector3d>> positions = readVector3dTupleFile(result_path+"/posChain_"+testName+"_"+std::to_string(ns)+".dat");
  cout<<"# of forces= "<<forces.size()<<endl;

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  // render window interacter
  // create renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(1., 1., 1.); // black
  renderWindow->AddRenderer(renderer);


  double maxForce = 0.0;
  double meanForce = 0.0;
  for(Vector3d force : forces) {
    maxForce = max(maxForce, force.norm());
    meanForce += force.norm();
  }
  meanForce /= forces.size();
  cout<<"at station "<<ns<<" max force is: "<<maxForce<<endl;

  //normalize such that the largest force has radius 5
  for(int i = 0; i < forces.size(); i+=1 ){
    if(i % 200 == 0){ cout<<"finishing "<<i<<" contact."<<endl; }
    auto [pos1, pos2] = positions[i];
    Vector3d position = (pos1+pos2)/2.;
    if(abs(pos1(0)) > 500. || abs(pos1(1)) > 500. || pos1(2) < 0. || pos1(2) > 1600. ){ continue; }
    //if(pos1(0) > 0.) continue;
    Vector3d force    = forces[i];
    Vector3d posDiff = pos1 - pos2;
    //cylinder
    vtkNew<vtkCylinderSource> cylinder;
    cylinder->SetHeight(posDiff.norm());                           // height
    cylinder->SetRadius(1.+2.*pow(force.norm()/maxForce, 1.));                            // radius
    cylinder->SetResolution(20);                        // number of slides
    cylinder->SetCenter(position(0), position(1), position(2));  // set center
    vtkNew<vtkPolyDataMapper> cylinderMapper;
    cylinderMapper->SetInputConnection(cylinder->GetOutputPort());
    vtkNew<vtkActor> clyinderActor;
    clyinderActor->SetMapper(cylinderMapper);
    if(force.norm() > 2.*meanForce){
      clyinderActor->GetProperty()->SetColor(1.,1-pow(force.norm()/maxForce, .25), 0.);
      clyinderActor->GetProperty()->SetOpacity(pow(force.norm()/maxForce,.5));
    }else{
      clyinderActor->GetProperty()->SetColor(0.,0.5, 1.);
      //clyinderActor->GetProperty()->SetColor(0.,pow(force.norm()/maxForce, 1.), 1.);
      clyinderActor->GetProperty()->SetOpacity(pow(force.norm()/maxForce,1.));
    }

    //clyinderActor->GetProperty()->SetColor(0.2, 0.2, 0.2);

    //clyinderActor->GetProperty()->SetOpacity(1.);
    clyinderActor->SetOrigin(position(0), position(1), position(2));
    double alpha = asin(posDiff(2)/posDiff.norm());
    double beta  = atan(posDiff(0)/posDiff(1));
    //cout<<"beta = "<<beta*180./3.1415926<<" alpha = "<<alpha*180.0/3.1415926<<endl;
    clyinderActor->RotateZ(-beta*180./3.1415926);
    if((posDiff(0) < 0) && (posDiff(1) < 0) || (posDiff(0) > 0) && (posDiff(1) < 0)){ alpha = -alpha; }
    clyinderActor->RotateX(alpha*180./3.1415926);

    renderer->AddActor(clyinderActor);
  }
  // Camera
  vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
  renderer->SetActiveCamera(camera);
  //Full View
  camera->SetViewUp(0,0,1);
  camera->SetPosition(300,0,300);
  camera->SetFocalPoint(0,0,250);
  camera->SetClippingRange(3,100);
  // THIS FUCKING THING RIGHT HERE (mak sure that the object is between
  //the first value and the second value wrt the camera's position/focal point)
  // THIS FUCKING THING RIGHT HERE (mak sure that the object is between
  //the first value and the second value wrt the camera's position/focal point)
  camera->Zoom(.9);

  // Render and interact
  renderWindow->Render();
  renderer->ResetCamera();

  vtkSmartPointer<vtkRenderLargeImage> renderLarge = vtkSmartPointer<vtkRenderLargeImage>::New();
  renderLarge->SetInput(renderer);
  renderLarge->SetMagnification(12);

  stringstream fname;
  fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/forceChain_"<<ns<<".png";

  //writer for saving later
  vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(fname.str().c_str());
	writer->SetInputConnection(renderLarge->GetOutputPort());
	writer->Write();

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Start();

	cout << "program terminated" << endl;
	return EXIT_SUCCESS;
}
