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
#include <vtkSphereSource.h>
#include <vtkPlaneSource.h>
#include <random>
#include "omp.h"
#include <algorithm>

int main(int argc, char *argv[])
{
  size_t ns = atoi(argv[1]);
	size_t ng = 42684; // number of grains 19330
  size_t total_ns = 40;
  string testName = "parallel_benchmark";
  string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/results";

  vector<Vector3d> allPositions = readPositionFile(result_path+"/positions_"+testName+".dat", ng*total_ns);
  vector<Vector4d> allRotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", ng*total_ns);

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  // render window interacter
  // create renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(1.0, 1.0, 1.0); // white
  renderWindow->AddRenderer(renderer);

  std::mt19937 gen(123123);
  std::uniform_real_distribution<double> rand_real(0., 1.0);

  vector<Vector3d> prevSteps = readPositionFile("/home/hasitha/Desktop/data/fabric/"+testName+"/first_4000_steps"+"/positions_"+testName+".dat", ng*total_ns);
  vector<Vector3d> initPositions(prevSteps.begin(), prevSteps.begin()+ng);
  vector<int> rightGrains;
  for(size_t i = 0; i < initPositions.size(); ++i){
    if(initPositions[i](1) > 300.) rightGrains.push_back(i);
  }

  vector<Vector3d> positions(allPositions.begin()+ns*ng, allPositions.begin()+(ns+1)*ng);
  vector<Vector4d> rotations(allRotations.begin()+ns*ng, allRotations.begin()+(ns+1)*ng);
  vector<int> grainList = {30173, 28612, 29064};

  for (size_t g = 0; g < ng; g+=1) {
    if(g % 200 == 0) cout << g << endl;
    Vector3d pos = positions[g];
    //if(std::find(grainList.begin(), grainList.end(), g) == grainList.end()) continue;
    stringstream fname;
    fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<g+1<<".dat";
    Polyhedron poly = readPolyFile(fname.str());
    vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, positions[g], rotations[g]/rotations[g].norm());
    // Create mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(trianglePolyData);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    if(std::find(rightGrains.begin(), rightGrains.end(), g) != rightGrains.end()){
      actor->GetProperty()->SetColor(1.0, 0., 0.);
    }else{
      actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
    }
    actor->GetProperty()->SetOpacity(1.);
    renderer->AddActor(actor);
  }

/*
{
  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetPoint1(400.,0.,0.);
  planeSource->SetPoint2(0.,400.,0.);
  planeSource->SetCenter(200.,200.,0.);
  planeSource->Update();
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(planeSource->GetOutputPort());
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(0.15, 0., 0.);
  actor->GetProperty()->SetOpacity(0.25);
  renderer->AddActor(actor);
}

{
  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetPoint1(400.,0.,0.);
  planeSource->SetPoint2(0.,0.,400.);
  planeSource->SetCenter(200.,400.,200.);
  planeSource->Update();
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(planeSource->GetOutputPort());
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(0.25, 0.25, 0.25);
  actor->GetProperty()->SetOpacity(0.25);
  renderer->AddActor(actor);
}

{
  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetPoint1(0.,400.,0.);
  planeSource->SetPoint2(0.,0.,400.);
  planeSource->SetCenter(400.,200.,200.);
  planeSource->Update();
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(planeSource->GetOutputPort());
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(0.5, 0.5, 0.5);
  actor->GetProperty()->SetOpacity(0.25);
  renderer->AddActor(actor);
}
*/
  // Camera
  vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
  renderer->SetActiveCamera(camera);
  //Full View
  camera->SetViewUp(0,0,1);
  camera->SetPosition(600,600,600);
  //camera->SetPosition(0,0,1500);
  camera->SetFocalPoint(300,300,300);
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
  fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/step_"<<ns<<".png";

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
