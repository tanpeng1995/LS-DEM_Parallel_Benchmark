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

int main(int argc, char *argv[])
{
  size_t total_ns = atoi(argv[1]);
  string testName = "7-6B_400";
  string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/results";

  double maxForce = 0.0;
  for (int ns = 0; ns < total_ns; ++ns){
    std::cout<<"currently processing "<<ns<<"-th step"<<std::endl;
    vector<Vector3d> forces    = readVectorFile(result_path+"/forceChain_"+testName+"_"+std::to_string(ns)+".dat");
  	vector<Vector3d> positions = readVectorFile(result_path+"/posChain_"+testName+"_"+std::to_string(ns)+".dat");

    for(Vector3d force : forces) {
      if(force.norm() > maxForce) {
        maxForce = force.norm();
      }
    }
  }





  for (int ns = 0; ns < total_ns; ++ns){
    std::cout<<"currently processing "<<ns<<"-th step"<<std::endl;
    vector<Vector3d> forces    = readVectorFile(result_path+"/forceChain_"+testName+"_"+std::to_string(ns)+".dat");
  	vector<Vector3d> positions = readVectorFile(result_path+"/posChain_"+testName+"_"+std::to_string(ns)+".dat");
    cout<<"# of forces= "<<forces.size()<<endl;

    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    // render window interacter
    // create renderer
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0., 0., 0.); // black
    renderWindow->AddRenderer(renderer);

    /*
    double maxForce = 0.0;
    for(Vector3d force : forces) {
      if(force.norm() > maxForce) {
        maxForce = force.norm();
      }
    }*/

    //normalize such that the largest force has radius 5
    for(int i = 0; i < forces.size(); ++i){

      Vector3d position = positions[i];
      Vector3d force    = forces[i];
      //cylinder
      vtkNew<vtkCylinderSource> cylinder;
      cylinder->SetHeight(20.0+40.0*(force.norm()/maxForce-0.5));                           // height
      cylinder->SetRadius(2.5*force.norm()/maxForce);                            // radius
      cylinder->SetResolution(100);                        // number of slides
      cylinder->SetCenter(position(0), position(1), position(2));  // set center
      vtkNew<vtkPolyDataMapper> cylinderMapper;
      cylinderMapper->SetInputConnection(cylinder->GetOutputPort());
      vtkNew<vtkActor> clyinderActor;
      clyinderActor->SetMapper(cylinderMapper);
      clyinderActor->GetProperty()->SetColor(1.0,1-(force.norm()/maxForce), 0.);
      clyinderActor->GetProperty()->SetOpacity(force.norm()/maxForce);
      //clyinderActor->GetProperty()->SetOpacity(1.);
      clyinderActor->SetOrigin(position(0), position(1), position(2));
      double alpha = asin(force(2)/force.norm());
      double beta  = atan(force(0)/force(1));
      //cout<<"beta = "<<beta*180./3.1415926<<" alpha = "<<alpha*180.0/3.1415926<<endl;
      clyinderActor->RotateZ(-beta*180./3.1415926);
      if((force(0) < 0) && (force(1) < 0) || (force(0) > 0) && (force(1) < 0)){ alpha = -alpha; }
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
    stringstream fname;
    fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/gforceChain_"<<ns<<".png";

    //writer for saving later
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
  	writer->SetFileName(fname.str().c_str());
  	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  	writer->Write();

    //vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    //renderWindowInteractor->SetRenderWindow(renderWindow);
    //renderWindowInteractor->Start();
  }

	cout << "program terminated" << endl;
	return EXIT_SUCCESS;
}
