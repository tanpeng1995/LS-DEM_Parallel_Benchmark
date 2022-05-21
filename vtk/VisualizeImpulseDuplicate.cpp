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
	size_t ng = 200;
  size_t ns = 80;
	size_t N = 2;
	size_t NG = ng * N * N * N;
	string testName = "uniform_impulse";
	string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/results";
	vector<Vector3d> allPositions = readPositionFile(result_path+"/positions_"+testName+".dat", NG*ns);
	vector<Vector4d> allRotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", NG*ns);

  for(size_t i = 0; i < ns; i+=4){
		cout<<"This is "<<i<<"-th frames"<<endl<<endl;
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    // render window interacter
    // create renderer
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1.0, 1.0, 1.0); // white
    renderWindow->AddRenderer(renderer);

    std::mt19937 gen(123123);
  	std::uniform_real_distribution<double> rand_real(0., 1.0);

    vector<Vector3d> positions(allPositions.begin()+i*NG, allPositions.begin()+(i+1)*NG);
    vector<Vector4d> rotations(allRotations.begin()+i*NG, allRotations.begin()+(i+1)*NG);

		for(int nz = 0; nz < N; nz++){
			for(int ny = 0; ny < N; ny++){
				for(int nx = 0; nx < N; nx++){
					int disp_i = (N*N*nz + N*ny + nx) * ng;
					for (size_t g = 0; g < ng; g+=1) {
						int gid = g + disp_i;
						if(gid % 20 == 0) cout<<gid<<endl;
						stringstream fname;
						double mass = 0;
				    fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<g+1<<".dat";
				    Polyhedron poly = readPolyFile2(fname.str());
			      vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, positions[gid], rotations[gid]/rotations[gid].norm());
			      // Create mapper and actor
			      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			      mapper->SetInputData(trianglePolyData);
			      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			      actor->SetMapper(mapper);
						actor->GetProperty()->SetColor(0.9,0.9,0.9);
						actor->GetProperty()->SetOpacity(1);
						renderer->AddActor(actor);
			    }
				}
			}
		}


		{
			vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
			planeSource->SetPoint1(800.,0.,0.);
			planeSource->SetPoint2(0.,800.,0.);
			planeSource->SetCenter(400.,400.,0.);
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
			planeSource->SetPoint1(800.,0.,0.);
			planeSource->SetPoint2(0.,0.,800.);
			planeSource->SetCenter(400.,0.,400.);
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
			planeSource->SetPoint1(0.,800.,0.);
			planeSource->SetPoint2(0.,0.,800.);
			planeSource->SetCenter(0.,400.,400.);
			planeSource->Update();
			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputConnection(planeSource->GetOutputPort());
			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(0.5, 0.5, 0.5);
			actor->GetProperty()->SetOpacity(0.25);
			renderer->AddActor(actor);
		}

    // Camera
  	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
  	renderer->SetActiveCamera(camera);

		camera->SetViewUp(0,0,1);
	  camera->SetPosition(100,100,100);
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
  	stringstream fname; fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/continue_4_"<<i<<".png";

		//writer for saving later
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
  	writer->SetFileName(fname.str().c_str());
  	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  	writer->Write();
  }

	cout << "program terminated" << endl;
	return EXIT_SUCCESS;
}
