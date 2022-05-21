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
	size_t ng = 71692;
  size_t ns = 130;
	string testName = "7-6B_cylinder";
	string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/results";
	vector<Vector3d> allPositions = readPositionFile(result_path+"/positions_"+testName+".dat", ng*ns);
	vector<Vector4d> allRotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", ng*ns);

	double ball_radius = 0;
	size_t nb = 0;
  vector<Vector3d> allBallPositions = readSphereMembrane(result_path+"/walls_"+testName+".dat", nb, ball_radius, ns);
	std::cout<<"ball_radius = "<<ball_radius<<std::endl;
  std::cout<<"#balls = "<<nb<<std::endl;

	double ball_radius2 = 0;
	size_t nb2 = 0;
  vector<Vector3d> allBallPositions2 = readSphereMembrane(result_path+"/external_walls_"+testName+".dat", nb2, ball_radius2, ns);
  std::cout<<"external ball_radius = "<<ball_radius2<<std::endl;
  std::cout<<"#externl balls = "<<nb2<<std::endl;

	std::vector<Polyhedron> grainPolys(ng);
	for (size_t g = 0; g < ng; g+=1) {
		if(g % 100 == 0){
			cout<<"finish constructing "<<g<<" polys."<<endl;
		}
		stringstream fname;
		fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<g+1<<".dat";
		grainPolys[g] = readPolyFile2(fname.str());
	}

  for(size_t i = 10; i < ns; i+=2){
		cout<<"This is "<<i<<"-th frames"<<endl<<endl;
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    // render window interacter
    // create renderer
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1.0, 1.0, 1.0); // white
    renderWindow->AddRenderer(renderer);

    std::mt19937 gen(123123);
  	std::uniform_real_distribution<double> rand_real(0., 1.0);
    vector<Vector3d> positions(allPositions.begin()+i*ng, allPositions.begin()+(i+1)*ng);
    vector<Vector4d> rotations(allRotations.begin()+i*ng, allRotations.begin()+(i+1)*ng);
		vector<size_t> badList = {182, 3819, 3820, 3901, 4096,4682, 5430, 5789, 5791, 5842, 5946, 6118, 6690,
			 6753, 6757, 6808, 7100, 7412, 7670, 9298, 25366, 32982, 34843, 65145, 65403, 65818, 66434, 66806,
			 67170, 67225, 67526, 68085, 69374, 69713, 69926, 70312, 70717, 70827, 71045, 71653};

    for (size_t g = 0; g < ng; g+=1) {
			if(std::find(badList.begin(), badList.end(), g) != badList.end()) continue;
			if(g % 100 == 0){
				cout<<"finish triangulating "<<g<<" polys. in station "<<i<<endl;
			}
      vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(grainPolys[g], positions[g], rotations[g]/rotations[g].norm());
      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputData(trianglePolyData);

      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(0.9,0.9,0.9);
			actor->GetProperty()->SetOpacity(1);
			renderer->AddActor(actor);
    }

/*
		{
		  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
		  planeSource->SetPoint1(700.,0.,0.);
		  planeSource->SetPoint2(0.,700.,0.);
		  planeSource->SetCenter(0.,0.,0.);
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
    vector<Vector3d> ballPositions(allBallPositions.begin()+i*nb, allBallPositions.begin()+(i+1)*nb);
    for(size_t b = 0; b < nb; b++){
      vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
      sphereSource->SetCenter(ballPositions[b](0), ballPositions[b](1), ballPositions[b](2));
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

  vector<Vector3d> ballPositions2(allBallPositions2.begin()+i*nb2, allBallPositions2.begin()+(i+1)*nb2);
/*
  for(size_t b = 0; b < nb2; b++){
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetCenter(ballPositions2[b](0), ballPositions2[b](1), ballPositions2[b](2));
    sphereSource->SetRadius(ball_radius2);
    sphereSource->Update();
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphereSource->GetOutputPort());
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
    actor->GetProperty()->SetOpacity(0.5);
		renderer->AddActor(actor);
  }
*/
    // Camera
  	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
  	renderer->SetActiveCamera(camera);

		camera->SetViewUp(0,0,1);
	  camera->SetPosition(300,-200,300);
	  camera->SetFocalPoint(0,0,250);
	  camera->SetClippingRange(3,100);
	  // THIS FUCKING THING RIGHT HERE (make sure that the object is between
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
  	stringstream fname; fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/good_shear_band_"<<i<<".png";

		//writer for saving later
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
  	writer->SetFileName(fname.str().c_str());
  	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  	writer->Write();
  }

	cout << "program terminated" << endl;
	return EXIT_SUCCESS;
}
