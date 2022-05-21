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
  size_t ns = 80;
	string testName = "impulse_400";
	string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/results";
	vector<Vector3d> allPositions = readPositionFile(result_path+"/positions_"+testName+".dat", ng*ns);
	vector<Vector4d> allRotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", ng*ns);

	double ball_radius = 0;
	size_t nb = 0;
  vector<Vector3d> allBallPositions = readSphereMembrane(result_path+"/walls_"+testName+".dat", nb, ball_radius, ns);
	std::cout<<"ball_radius = "<<ball_radius<<std::endl;
  std::cout<<"#balls = "<<nb<<std::endl;

  for(size_t i = 0; i < ns; i+=1){
		cout<<"This is "<<i<<"-th frames"<<endl<<endl;
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    // render window interacter
    // create renderer
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1.0, 1.0, 1.0); // white
    renderWindow->AddRenderer(renderer);

		for(size_t b = 0; b < nb; b+=1){
	    if(b % 100 == 0) cout<<b<<endl;
	    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
	    int pos = i*nb+b;
			Vector3d ballPos = allBallPositions[pos];
			if(ballPos(0)<0. || ballPos(0)>1800. || ballPos(1)<0. || ballPos(1)>600. || ballPos(2)<0.|| ballPos(2)>1800.) continue;
	    sphereSource->SetCenter(allBallPositions[pos](0), allBallPositions[pos](1), allBallPositions[pos](2));
	    //if(b % 1000 == 0) std::cout<<ballPositions[b]<<std::endl;
	    sphereSource->SetRadius(ball_radius);
	    sphereSource->Update();

	    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	    mapper->SetInputConnection(sphereSource->GetOutputPort());

	    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	    actor->SetMapper(mapper);
	    actor->GetProperty()->SetColor(1.0, 0., 0.);
	    actor->GetProperty()->SetOpacity(0.35);

	    renderer->AddActor(actor);
	  }


    std::mt19937 gen(123123);
  	std::uniform_real_distribution<double> rand_real(0., 1.0);

    vector<Vector3d> positions(allPositions.begin()+i*ng, allPositions.begin()+(i+1)*ng);
    vector<Vector4d> rotations(allRotations.begin()+i*ng, allRotations.begin()+(i+1)*ng);

		//vector<int> grainList = {30, 99};

	  for (size_t g = 0; g < ng; g+=1) {
	    Vector3d pos = positions[g];
			if(pos(0)<0. || pos(0)>1800. || pos(1)<0. || pos(1)>600. || pos(2)<0.|| pos(2)>1800.) continue;
	    //if(std::find(grainList.begin(), grainList.end(), g) == grainList.end() ) continue;
	    if(g % 20 == 0) cout<<g<<endl;
	    stringstream fname;
	    double mass = 0;
	    fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<g+1<<".dat";
	    Polyhedron poly = readPolyFile2(fname.str());
	    vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, positions[g], rotations[g]/rotations[g].norm());
	    // Create mapper and actor
	    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	    mapper->SetInputData(trianglePolyData);

	    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	    actor->SetMapper(mapper);
	    actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
	    actor->GetProperty()->SetOpacity(1.);
	    renderer->AddActor(actor);
	  }

		{
		  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
			planeSource->SetPoint1(1035., 0., -598.);
		  planeSource->SetPoint2(0.,600.,0.);
		  planeSource->SetCenter(990., 300., 468.);
		  planeSource->Update();
		  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		  mapper->SetInputConnection(planeSource->GetOutputPort());
		  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		  actor->SetMapper(mapper);
		  actor->GetProperty()->SetColor(0.15, 0., 0.);
		  actor->GetProperty()->SetOpacity(0.25);
		  renderer->AddActor(actor);
		}

    // Camera
  	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
  	renderer->SetActiveCamera(camera);

		camera->SetViewUp(0,0,1);
	  camera->SetPosition(900,-300,300);
	  camera->SetFocalPoint(900,300,250);
	  camera->SetClippingRange(3,100);
	  // THIS FUCKING THING RIGHT HERE (mak sure that the object is between
	  //the first value and the second value wrt the camera's position/focal point)
	  camera->Zoom(0.9);

  	// Render and interact
  	renderWindow->Render();
  	renderer->ResetCamera();

  	//Screenshot
  	vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
  	windowToImageFilter->SetInput(renderWindow);
  	windowToImageFilter->SetScale(3,3);
  	windowToImageFilter->ReadFrontBufferOff();
  	windowToImageFilter->Update();
  	stringstream fname; fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/paper_picture_"<<i<<".png";

		//writer for saving later
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
  	writer->SetFileName(fname.str().c_str());
  	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  	writer->Write();
  }

	cout << "program terminated" << endl;
	return EXIT_SUCCESS;
}
