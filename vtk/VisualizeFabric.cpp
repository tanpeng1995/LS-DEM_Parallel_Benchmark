/*
 * VisualizeFabric.cpp
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
#include <vtkRenderLargeImage.h>

int main(int argc, char *argv[])
{
	size_t ng = 19606;
	size_t gg = atoi(argv[1]);
	string testName = "cylinder_belgium";
	string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/InitState";
	vector<Vector3d> positions = readPositionFile(result_path+"/positions_"+testName+".dat", ng);
	vector<Vector4d> rotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", ng);

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  // render window interacter
  // create renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(1.0, 1.0, 1.0); // white
  renderWindow->AddRenderer(renderer);

  std::mt19937 gen(123123);
	std::uniform_real_distribution<double> rand_real(0., 1.0);

  for (size_t g = 0; g < ng; g+=1) {
    if(g % 10 == 0) std::cout<<g<<std::endl;
		//if(g != gg) continue;
		Vector3d pos = positions[g];
		stringstream fname;
    fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<g+1<<".dat";
		double mass;
    Polyhedron poly = readPolyFile2(fname.str());
    vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, pos, rotations[g]/rotations[g].norm());
    // Create mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(trianglePolyData);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    renderer->AddActor(actor);
    //actor->GetProperty()->SetColor(rand_real(gen), rand_real(gen), rand_real(gen));
		actor->GetProperty()->SetColor(0.9,0.9,0.9);
    actor->GetProperty()->SetOpacity(1.);
  }

/*
// BOTTOM
	{
		vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
		planeSource->SetPoint1(600.,0.,0.);
		planeSource->SetPoint2(0.,600.,0.);
		planeSource->SetCenter(300.,300.,0.);
		planeSource->Update();
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(planeSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
		actor->GetProperty()->SetOpacity(0.15);
		renderer->AddActor(actor);
	}

//UP
{
	vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
	planeSource->SetPoint1(600.,0.,0.);
	planeSource->SetPoint2(0.,600.,0.);
	planeSource->SetCenter(300.,300.,600.);
	planeSource->Update();
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(planeSource->GetOutputPort());
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
	actor->GetProperty()->SetOpacity(0.15);
	renderer->AddActor(actor);
}

//UP
{
	vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
	planeSource->SetPoint1(600.,0.,0.);
	planeSource->SetPoint2(0.,600.,0.);
	planeSource->SetCenter(300.,300.,300.);
	planeSource->Update();
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(planeSource->GetOutputPort());
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(0.5, 0.5, 0.5);
	actor->GetProperty()->SetOpacity(0.15);
	renderer->AddActor(actor);
}

//LEFT
	{
		vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
		planeSource->SetPoint1(600.,0.,0.);
		planeSource->SetPoint2(0.,0.,600.);
		planeSource->SetCenter(300.,600.,300.);
		planeSource->Update();
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(planeSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
		actor->GetProperty()->SetOpacity(0.15);
		renderer->AddActor(actor);
	}

//RIGHT
	{
		vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
		planeSource->SetPoint1(600.,0.,0.);
		planeSource->SetPoint2(0.,0.,600.);
		planeSource->SetCenter(300.,0.,300.);
		planeSource->Update();
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(planeSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
		actor->GetProperty()->SetOpacity(0.15);
		renderer->AddActor(actor);
	}

	{
		vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
		planeSource->SetPoint1(600.,0.,0.);
		planeSource->SetPoint2(0.,0.,600.);
		planeSource->SetCenter(300.,300.,300.);
		planeSource->Update();
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(planeSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.5, 0.5, 0.5);
		actor->GetProperty()->SetOpacity(0.15);
		renderer->AddActor(actor);
	}


//BACK
	{
		vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
		planeSource->SetPoint1(0.,600.,0.);
		planeSource->SetPoint2(0.,0.,600.);
		planeSource->SetCenter(600.,300.,300.);
		planeSource->Update();
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(planeSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
		actor->GetProperty()->SetOpacity(0.15);
		renderer->AddActor(actor);
	}
//FRONT
{
	vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
	planeSource->SetPoint1(0.,600.,0.);
	planeSource->SetPoint2(0.,0.,600.);
	planeSource->SetCenter(0.,300.,300.);
	planeSource->Update();
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(planeSource->GetOutputPort());
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
	actor->GetProperty()->SetOpacity(0.15);
	renderer->AddActor(actor);
}

{
	vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
	planeSource->SetPoint1(0.,600.,0.);
	planeSource->SetPoint2(0.,0.,600.);
	planeSource->SetCenter(300.,300.,300.);
	planeSource->Update();
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(planeSource->GetOutputPort());
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(0.5, 0.5, 0.5);
	actor->GetProperty()->SetOpacity(0.15);
	renderer->AddActor(actor);
}
*/
  // Camera
	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
	renderer->SetActiveCamera(camera);
	//Full View
	camera->SetViewUp(0,0,1);
	camera->SetPosition(0,0,150);
	//camera->SetPosition(0,0,1500);
	camera->SetFocalPoint(300,300,0);
	camera->SetClippingRange(3,100);
	// THIS FUCKING THING RIGHT HERE (mak sure that the object is between
	//the first value and the second value wrt the camera's position/focal point)
	camera->Zoom(1.1);

	// Render and interact
	renderWindow->Render();
	renderer->ResetCamera();

	vtkSmartPointer<vtkRenderLargeImage> renderLarge = vtkSmartPointer<vtkRenderLargeImage>::New();
  renderLarge->SetInput(renderer);
  renderLarge->SetMagnification(2);

	stringstream fname;
	fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/Fabric_new3.png";

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
