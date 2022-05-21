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

vector<int> readLargeRotationFile(string filename){
  ifstream file(filename.c_str());
  string line;
  vector<int> integers;
  while (getline(file, line)) {
		integers.push_back( atoi(line.c_str()) );
	}
  cout<<"large rotation grains #: "<<integers.size()<<endl;
	return integers;
}

vector<double> readRotationAngles(string filename, double & maxRot, double & meanRot){
  ifstream file(filename.c_str());
  string line;
  vector<double> angles;
  maxRot = 0;
  meanRot = 0.;
  while (getline(file, line)) {
    double rot = abs( atof(line.c_str()) );
    maxRot = max( maxRot, rot );
    meanRot += rot;
		angles.push_back( rot );
	}
  meanRot /= angles.size();
  cout<<"rotation angles #: "<<angles.size()<<" mean: "<<meanRot<<" max: "<<maxRot<<endl;
	return angles;
}

int main(int argc, char *argv[])
{
  size_t ns = atoi(argv[1]);
	size_t ng = 19606;
  size_t total_ns = 200;
  string testName = "cylinder_belgium";
  string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/results";

	vector<Vector3d> allPositions = readPositionFile(result_path+"/positions_"+testName+".dat", ng*total_ns);
	vector<Vector4d> allRotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", ng*total_ns);
  cout<<"allPositions size = "<<allPositions.size()<<endl;
  double ball_radius = 0;
	size_t nb = 0;
  vector<Vector3d> allBallPositions = readSphereMembrane(result_path+"/walls_"+testName+".dat", nb, ball_radius, total_ns);
  std::cout<<"ball_radius = "<<ball_radius<<std::endl;
  std::cout<<"#balls = "<<nb<<std::endl;

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(1.0, 1.0, 1.0); // white
  renderWindow->AddRenderer(renderer);
/*
  for(size_t b = 0; b < nb; b+=1){
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    int pos = ns*nb+b;
    Vector3d ballPos = allBallPositions[pos];
    //if(ballPos(0) > 0.) continue;
    sphereSource->SetCenter(ballPos(0), ballPos(1), ballPos(2));
    sphereSource->SetRadius(ball_radius);
    sphereSource->Update();
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphereSource->GetOutputPort());
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
    actor->GetProperty()->SetOpacity(0.25);
    renderer->AddActor(actor);
  }
*/


  double ball_radius2 = 0;
	size_t nb2 = 0;
  vector<Vector3d> allBallPositions2 = readSphereMembrane(result_path+"/external_walls_"+testName+".dat", nb2, ball_radius2, total_ns);
  std::cout<<"external ball_radius = "<<ball_radius2<<std::endl;
  std::cout<<"#externl balls = "<<nb2<<std::endl;
/*
  for(size_t b = 0; b < nb2; b++){
    if(b % 100 == 0) cout<<b<<endl;
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
    actor->GetProperty()->SetOpacity(0.075);

    renderer->AddActor(actor);
  }
*/

  std::mt19937 gen(123123);
  std::uniform_real_distribution<double> rand_real(-1.0, 1.0);

  vector<Vector3d> positions(allPositions.begin()+ns*ng, allPositions.begin()+(ns+1)*ng);
  vector<Vector4d> rotations(allRotations.begin()+ns*ng, allRotations.begin()+(ns+1)*ng);
  vector<int> largeRotationGrain = readLargeRotationFile(result_path+"/large_rotation.dat");

  double maxRot, meanRot;
  vector<double> rotationAngles = readRotationAngles(result_path+"/rotation_angles.dat", maxRot, meanRot);

  double largem = 0;
  for(int i : largeRotationGrain){
    largem += rotationAngles[i];
  }
  largem /= largeRotationGrain.size();


  vector<int> badList;
  Vector3d w1;
  w1 << 0.335, -0.31,  -0.89;
  double w[3] = { 0.335, -0.31,  -0.89 };
  double height = 600.;
  double dplane = -w1(2)*height;

  for (size_t g = 0; g < ng; g+=1) {
    Vector3d pos = positions[g];
    Vector4d rot = rotations[g];
    bool flag = true;
    if(g % 100 == 0) cout<<g<<endl;
    //if( abs(pos.dot(w1)+height) > 300.) continue;
    if(pos(0) > 0.) continue;
    //if(sqrt(pos(0)*pos(0) + pos(1)*pos(1)) < 140.) continue;
    if(abs(pos(0)) > 600. || abs(pos(1)) > 600. || pos(2) < 0. || pos(2) > 1700. ){ badList.push_back(g); continue; }
    if(std::find(largeRotationGrain.begin(), largeRotationGrain.end(), g) == largeRotationGrain.end() ){ flag = false;}
    //if(g >= 8000 && g <= 11000 && g % 4 != 0) { rotationAngles[g] /= 2; flag = false; continue; }
    stringstream fname;
    fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<g+1<<".dat";
    Polyhedron poly = readPolyFile2(fname.str());
    vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, positions[g], rotations[g]/rotations[g].norm());
    // Create mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(trianglePolyData);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    if(flag){
      actor->GetProperty()->SetColor(1., pow(1.-(rotationAngles[g]/100.), .25), 0.);
      actor->GetProperty()->SetOpacity(pow((rotationAngles[g]/100.), 2.));
      //actor->GetProperty()->SetOpacity(1.);
    }else{
      actor->GetProperty()->SetColor(0., pow(rotationAngles[g]/meanRot/1.5, 2.), 1. );
      actor->GetProperty()->SetOpacity(.75);
    }
    renderer->AddActor(actor);
  }

  for(int g : badList){
    cout<<g<<" ";
  }
  cout<<endl;
/*
{
  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetPoint1(800.,0.,0.);
  planeSource->SetPoint2(0.,800.,0.);
  planeSource->SetCenter(0.,0.,height);
  planeSource->SetNormal(w);
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
/*
{
  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetPoint1(400.,0.,0.);
  planeSource->SetPoint2(0.,400.,0.);
  planeSource->SetCenter(0.,0.,0.);
  planeSource->Update();
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(planeSource->GetOutputPort());
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(0.3, 0.3, 0.3);
  actor->GetProperty()->SetOpacity(0.2);
  renderer->AddActor(actor);
}
*/

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

  vtkSmartPointer<vtkRenderLargeImage> renderLarge = vtkSmartPointer<vtkRenderLargeImage>::New();
  renderLarge->SetInput(renderer);
  renderLarge->SetMagnification(12);

  stringstream fname;
  fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/test_"<<ns<<".png";

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
