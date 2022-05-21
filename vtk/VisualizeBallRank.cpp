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
#include <vtkRenderLargeImage.h>

int main(int argc, char *argv[])
{
  size_t ns = atoi(argv[1]);
  size_t rank = atoi(argv[2]);
  size_t total_ns = 150;
  string testName = "7-6B_cylinder";
  string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/results";

	vector<std::set<int>> ballRanks = readSetFile(result_path+"/ballRank_"+testName+"_"+std::to_string(ns)+".dat", rank);

  double ball_radius = 0;
	size_t nb = 0;
  vector<Vector3d> allBallPositions = readSphereMembrane(result_path+"/walls_"+testName+".dat", nb, ball_radius, total_ns);
  std::cout<<"ball_radius = "<<ball_radius<<std::endl;
  std::cout<<"#balls = "<<nb<<std::endl;

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(1.0, 1.0, 1.0); // white
  renderWindow->AddRenderer(renderer);

  for(size_t b = 0; b < nb; b+=1){
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    int pos = ns*nb+b;
    Vector3d ballPos = allBallPositions[pos];
    sphereSource->SetCenter(allBallPositions[pos](0), allBallPositions[pos](1), allBallPositions[pos](2));
    sphereSource->SetRadius(ball_radius);
    sphereSource->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphereSource->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
    actor->GetProperty()->SetOpacity(0.1);
    renderer->AddActor(actor);
  }


  std::mt19937 gen(123123);
  std::uniform_real_distribution<double> rand_real(0., 1.0);
  for(size_t i = 0; i < rank; ++i){
    double color[3] = {rand_real(gen), rand_real(gen), rand_real(gen)};
    cout<<"color for rank: "<<i<<" "<<color[0]<<" "<<color[1]<<" "<<color[2]<<endl;
    for(auto bid : ballRanks[i]){
      {
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        int pos = ns*nb+bid;
        Vector3d ballPos = allBallPositions[pos];
        sphereSource->SetCenter(allBallPositions[pos](0), allBallPositions[pos](1), allBallPositions[pos](2));
        sphereSource->SetRadius(ball_radius);
        sphereSource->Update();

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphereSource->GetOutputPort());

        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(color[0], color[1], color[2]);
        actor->GetProperty()->SetOpacity(0.95);
        renderer->AddActor(actor);
      }
    }
  }

  /*
  double ball_radius2 = 0;
	size_t nb2 = 0;
  vector<Vector3d> allBallPositions2 = readSphereMembrane(result_path+"/external_walls_"+testName+".dat", nb2, ball_radius2, total_ns);
  std::cout<<"external ball_radius = "<<ball_radius2<<std::endl;
  std::cout<<"#externl balls = "<<nb2<<std::endl;

  for(size_t b = 0; b < nb2; b++){
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    int pos = ns*nb2+b;
    sphereSource->SetCenter(allBallPositions2[pos](0), allBallPositions2[pos](1), allBallPositions2[pos](2));
    //if(b % 1000 == 0) std::cout<<ballPositions2[b]<<std::endl;
    sphereSource->SetRadius(ball_radius2);
    sphereSource->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphereSource->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.5, 0.5, 0.5);
    actor->GetProperty()->SetOpacity(0.05);

    renderer->AddActor(actor);
  }
  */

/* floor */
{
  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetPoint1(800.,0.,0.);
  planeSource->SetPoint2(0.,800.,0.);
  planeSource->SetCenter(0.,0.,0.);
  planeSource->Update();
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(planeSource->GetOutputPort());
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(0.25, 0.25, 0.25);
  actor->GetProperty()->SetOpacity(0.25);
  renderer->AddActor(actor);
}
/*
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
  camera->SetViewUp(0,0,1);
  camera->SetPosition(300,300,275);
  //camera->SetPosition(0,0,1500);
  camera->SetFocalPoint(0,0,250);
  camera->SetClippingRange(3,100);
  camera->Zoom(.9);

  // Render and interact
  renderWindow->Render();
  renderer->ResetCamera();

  //Screenshot
  /*
  vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowToImageFilter->SetInput(renderWindow);
  windowToImageFilter->SetScale(3,3);
  windowToImageFilter->ReadFrontBufferOff();
  windowToImageFilter->Update();
  */
  vtkSmartPointer<vtkRenderLargeImage> renderLarge = vtkSmartPointer<vtkRenderLargeImage>::New();
  renderLarge->SetInput(renderer);
  renderLarge->SetMagnification(8);

  stringstream fname;
  fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/ball_rank_"<<ns<<".png";

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
