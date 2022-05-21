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
#include "omp.h"

bool qualifiedGrain(const Vector3d & position, const Vector4d & quat,
	const string & morphfile, int gid, const double & lowHeight,
	const double & highHeight, const double & radius, const double & minMass){
	ifstream file(morphfile.c_str());
	if(file.fail()){
		cout << "Grain "<< gid + 1 << " does not exist." << endl;
		return false;
	}
	string line;
	string partial;
	istringstream iss;
	//mass
	getline(file, line);
	double mass = atof(line.c_str());
	if(mass < minMass) {
		return false;
	}
	// moment of inertia
	getline(file, line);
	//cmLset
	getline(file, line);
	// npoints
	getline(file, line);
	int npoints = atoi(line.c_str());
	// points
	Matrix3d R = Rfromq(quat);
	getline(file, line);
	Vector3d point;
	iss.str(line);
	for(int i = 0; i < npoints; ++i){
		getline(iss, partial, ' ');
		point(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(2) = atof(partial.c_str());
		point = R * point + position;
		if(point(2) < lowHeight || point(2) > highHeight ||
      point(0)*point(0) + point(1)*point(1) > radius*radius){
			return false;
		}
	}
	iss.clear();
	return true;
}

int main(int argc, char *argv[])
{
  size_t ns = atoi(argv[1]);
	size_t ng = 58377; // number of grains 19330
  size_t total_ns = 1;
  string testName = "7-6B_full";
  string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/InitState";

	vector<Vector3d> allPositions = readPositionFile(result_path+"/positions_"+testName+".dat", ng*total_ns);
  vector<Vector3d> initPostions(allPositions.begin(), allPositions.begin()+ng);
	vector<Vector4d> allRotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", ng*total_ns);
  vector<Vector4d> initRotations(allRotations.begin(), allRotations.begin()+ng);

  double ball_radius = 0;
	size_t nb = 0;
  vector<Vector3d> allBallPositions = readSphereMembrane(result_path+"/walls_"+testName+".dat", nb, ball_radius, total_ns);
  std::cout<<"ball_radius = "<<ball_radius<<std::endl;
  std::cout<<"#balls = "<<nb<<std::endl;

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  // render window interacter
  // create renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(1.0, 1.0, 1.0); // white
  renderWindow->AddRenderer(renderer);

  for(size_t b = 0; b < nb; b+=1){
    if(b % 100 == 0) cout<<b<<endl;
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    int pos = ns*nb+b;
    sphereSource->SetCenter(allBallPositions[pos](0), allBallPositions[pos](1), allBallPositions[pos](2));
    //if(b % 1000 == 0) std::cout<<ballPositions[b]<<std::endl;
    sphereSource->SetRadius(ball_radius);
    sphereSource->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphereSource->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
    actor->GetProperty()->SetOpacity(0.5);

    renderer->AddActor(actor);
  }


  double ball_radius2 = 0;
	size_t nb2 = 0;
  vector<Vector3d> allBallPositions2 = readSphereMembrane(result_path+"/external_walls_"+testName+".dat", nb2, ball_radius2, total_ns);
  std::cout<<"external ball_radius = "<<ball_radius2<<std::endl;
  std::cout<<"#externl balls = "<<nb2<<std::endl;
/*
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
    actor->GetProperty()->SetOpacity(0.25);

    renderer->AddActor(actor);
  }*/

  std::mt19937 gen(123123);
  std::uniform_real_distribution<double> rand_real(0., 1.0);

  vector<Vector3d> positions(allPositions.begin()+ns*ng, allPositions.begin()+(ns+1)*ng);
  vector<Vector4d> rotations(allRotations.begin()+ns*ng, allRotations.begin()+(ns+1)*ng);

  for (size_t g = 0; g < ng; g+=1) {
    if(g % 1000 == 0) cout<<g<<endl;

    //stringstream mname;
		//mname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Morphologies/morph_"<<g+1<<".dat";
    //if( !qualifiedGrain(positions[g], rotations[g], mname.str(), g, 0., 800., 200., 400.) ) continue;
    stringstream fname;
    fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<g+1<<".dat";
    ifstream file(fname.str().c_str());
    string   line;
    getline(file, line);
    double mass = atof(line.c_str());
    //cout<<"g = "<<g<<", mass = "<<mass<<endl;
    if(mass >= 400. ){
      Polyhedron poly = readPolyFile2(fname.str());
      vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, positions[g], rotations[g]/rotations[g].norm());
      // Create mapper and actor
      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputData(trianglePolyData);

      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);

      renderer->AddActor(actor);
      //actor->GetProperty()->SetColor(rand_real(gen), rand_real(gen), rand_real(gen));
      actor->GetProperty()->SetColor(0.9,0.9,0.9);
      actor->GetProperty()->SetOpacity(1);
    }else{

      vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
      sphereSource->SetCenter(positions[g](0), positions[g](1), positions[g](2));
      sphereSource->SetRadius(5.);
      sphereSource->Update();

      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputConnection(sphereSource->GetOutputPort());

      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);
      actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
      actor->GetProperty()->SetOpacity(1);

      renderer->AddActor(actor);

    }
  }




  // Camera
  vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
  renderer->SetActiveCamera(camera);
  //Full View
  camera->SetViewUp(0,0,1);
  camera->SetPosition(-1,1,100);
  //camera->SetPosition(0,0,1500);
  camera->SetFocalPoint(0,0,100);
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
  fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/one_station_"<<ns<<"_3.png";

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
