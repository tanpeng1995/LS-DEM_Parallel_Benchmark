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
#include <random>

int main(int argc, char *argv[])
{
  size_t ng = 3644;
  string testName = "fabric_600";
  string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/results";
  vector<Vector3d> contacts = readPositionFile(result_path+"/qhull_"+testName+".dat");
  vector<Vector3d> fullContacts = readPositionFile(result_path+"/flexible_"+testName+".dat");
  vector<Vector3d> positions = readPositionFile(result_path+"/positions_"+testName+".dat", ng);
	vector<Vector4d> rotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", ng);

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  // render window interacter
  // create renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(1.0, 1.0, 1.0); // white
  renderWindow->AddRenderer(renderer);
/*
  for (size_t g = 0; g < ng; g++) {
    if(g % 100 == 0){ std::cout<<g<<std::endl; }
		stringstream fname;
		fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<g+1<<".dat";
		Polyhedron poly = readPolyFile(fname.str());
		vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, positions[g], rotations[g]/rotations[g].norm());
		// Create mapper and actor
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(trianglePolyData);

		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		// thress largest particles are 8, 61, 66;
		actor->GetProperty()->SetColor(0.9,0.9,0.9);
		actor->GetProperty()->SetOpacity(0.15);

		renderer->AddActor(actor);
  }*/

  cout<<"#of quick hull is: "<<contacts.size()<<endl;
  for(size_t b = 0; b < contacts.size(); b+=1){
    if(b%500 == 0) cout<<b<<endl;
    //if(contacts[b](0) > 0. || contacts[b](2) < 100. ) continue;
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetCenter(contacts[b](0), contacts[b](1), contacts[b](2));
    //if(b % 1000 == 0) std::cout<<ballPositions[b]<<std::endl;
    sphereSource->SetRadius(2.5);
    sphereSource->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphereSource->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0., 0., 0.);
    actor->GetProperty()->SetOpacity(0.8);

    renderer->AddActor(actor);
  }
/*
  for(size_t b = 0; b < fullContacts.size(); b+=5){
    if(b%500 == 0) cout<<b<<endl;
    if(fullContacts[b](2) < 100. || fullContacts[b](2) > 1000.) continue;
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetCenter(fullContacts[b](0), fullContacts[b](1), fullContacts[b](2));
    //if(b % 1000 == 0) std::cout<<ballPositions[b]<<std::endl;
    sphereSource->SetRadius(2.5);
    sphereSource->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphereSource->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1.0, 0.0, 0.5);
    actor->GetProperty()->SetOpacity(0.8);

    renderer->AddActor(actor);
  }*/

  // Camera
  vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
  renderer->SetActiveCamera(camera);
  //Full View
  camera->SetViewUp(0,0,1);
  camera->SetPosition(-100,-30,0);
  //camera->SetPosition(0,0,1500);
  camera->SetFocalPoint(0,0,10);
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
  fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/flexibleContact.png";

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
