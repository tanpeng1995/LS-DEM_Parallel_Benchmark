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
#include <memory>

int main(int argc, char *argv[])
{
	size_t ng = 853;
  size_t total_ns = 120;
  string testName = "topography_long";
  string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/w_double_net";

	vector<Vector3d> allPositions = readPositionFile(result_path+"/positions_"+testName+".dat", ng*total_ns);
	vector<Vector4d> allRotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", ng*total_ns);

	double ball_radius = 0;
	size_t nb = 0;
	vector<Vector3d> allBallPositions = readSphereMembrane(result_path+"/walls_1_"+testName+".dat", nb, ball_radius, total_ns);
	std::cout<<"ball_radius for 1st wallballs = "<<ball_radius<<std::endl;
	std::cout<<"#balls = "<<nb<<std::endl;
	double ball_radius2 = 0;
	size_t nb2 = 0;
	vector<Vector3d> allBallPositions2 = readSphereMembrane(result_path+"/walls_2_"+testName+".dat", nb2, ball_radius2, total_ns);
	std::cout<<"ball_radius for 2nd wallballs = "<<ball_radius2<<std::endl;
	std::cout<<"#balls = "<<nb2<<std::endl;

  stringstream fname;
  fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/poly_topography_zoomin.dat";
  Polyhedron poly = readTopographyFile(fname.str());
  vtkSmartPointer<vtkPolyData> topographyPolyData = constructVtkPoly(poly, Vector3d(0.,0.,400.), 0.2247);

	double maxMass = std::numeric_limits<double>::min();
  std::vector<double> grainMass(ng);
  std::vector<Polyhedron> grainPolys(ng);
	for (size_t g = 0; g < ng; g+=1) {
		if(g % 20 == 0){
			cout<<"finish constructing "<<g<<" polys."<<endl;
		}
		stringstream fname; double mass;
		fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<g+1<<".dat";
		grainPolys[g] = readPolyFile3(fname.str(), mass);
    grainMass[g] = mass;
    maxMass = max(maxMass, mass);
	}

  for(int ns = 90; ns < total_ns; ns += 2){
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1.0, 1.0, 1.0); // white
    renderWindow->AddRenderer(renderer);

    /* digital topography */
    {
      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputData(topographyPolyData);

      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);
      actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
      actor->GetProperty()->SetOpacity(1.);
      renderer->AddActor(actor);
    }

		{
		  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
			planeSource->SetPoint1(1100.,0.,0.);
		  planeSource->SetPoint2(0.,800.,0.);
		  planeSource->SetCenter(1600.,320.,200.);
			//planeSource->SetNormal(0.,-0.052335956,0.9986295348);
		  planeSource->Update();
		  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		  mapper->SetInputConnection(planeSource->GetOutputPort());
		  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		  actor->SetMapper(mapper);
		  actor->GetProperty()->SetColor(1., 1., 1.);
		  actor->GetProperty()->SetOpacity(1.);
		  renderer->AddActor(actor);
		}

    cout<<"Constructing "<<ns<<"-th images"<<endl;

    vector<Vector3d> positions(allPositions.begin()+ns*ng, allPositions.begin()+(ns+1)*ng);
    vector<Vector4d> rotations(allRotations.begin()+ns*ng, allRotations.begin()+(ns+1)*ng);

    for (size_t g = 0; g < ng; g+=1) {
      vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(grainPolys[g], positions[g], rotations[g]/rotations[g].norm());
      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputData(trianglePolyData);

      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(1., 1. - grainMass[g]/maxMass, 0.);
      actor->GetProperty()->SetOpacity(1.);
      renderer->AddActor(actor);
    }

		for(size_t b = 0; b < nb; b+=1){
	    if(b % 100 == 0) cout<<b<<endl;
	    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
	    int pos = ns*nb+b;
	    Vector3d ballPos = allBallPositions[pos];
	    sphereSource->SetCenter(ballPos(0), ballPos(1), ballPos(2));
	    sphereSource->SetRadius(ball_radius);
	    sphereSource->Update();

	    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	    mapper->SetInputConnection(sphereSource->GetOutputPort());

	    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	    actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(0.307569, 0.814542, 0.454144);
	    actor->GetProperty()->SetOpacity(.5);
	    renderer->AddActor(actor);
	  }

		for(size_t b = 0; b < nb2; b+=1){
	    if(b % 100 == 0) cout<<b<<endl;
	    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
	    int pos = ns*nb2+b;
	    Vector3d ballPos = allBallPositions2[pos];
	    sphereSource->SetCenter(ballPos(0), ballPos(1), ballPos(2));
	    sphereSource->SetRadius(ball_radius2);
	    sphereSource->Update();

	    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	    mapper->SetInputConnection(sphereSource->GetOutputPort());

	    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	    actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(0.34266, 0.406228, 0.684375);
	    actor->GetProperty()->SetOpacity(.5);
	    renderer->AddActor(actor);
	  }


    // Camera
    vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
    renderer->SetActiveCamera(camera);
    //Full View
		camera->SetViewUp(0,0,1);
	  camera->SetPosition(1500,0,750);
	  camera->SetFocalPoint(750,320,150);
	  camera->SetClippingRange(3,100);
    // THIS FUCKING THING RIGHT HERE (mak sure that the object is between
    //the first value and the second value wrt the camera's position/focal point)
    camera->Zoom(.9);

    // Render and interact
    renderWindow->Render();
    renderer->ResetCamera();

    vtkSmartPointer<vtkRenderLargeImage> renderLarge = vtkSmartPointer<vtkRenderLargeImage>::New();
    renderLarge->SetInput(renderer);
    renderLarge->SetMagnification(4);

    stringstream fname;
    fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/topography_w_double_net_"<<ns<<".png";

    //writer for saving later
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
  	writer->SetFileName(fname.str().c_str());
  	writer->SetInputConnection(renderLarge->GetOutputPort());
  	writer->Write();
  }

	cout << "program terminated" << endl;
	return EXIT_SUCCESS;
}
