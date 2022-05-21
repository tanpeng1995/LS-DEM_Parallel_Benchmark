#include "definitions.h"
#include "readInputMethods.h"
#include "buildMethods.h"
#include <vtkPlaneSource.h>
#include  <vtkSphereSource.h>
#include <random>
#include "omp.h"

int main(int argc, char *argv[])
{
	size_t ns = atoi(argv[1]);
	size_t ng = 19377;
  size_t total_ns = 10;
	size_t N = 2;
	size_t NG = ng * N * N * N;
	string testName = "impulse_parallel";
	string result_path = "/home/hasitha/Desktop/data/fabric/"+testName+"/results";
	vector<Vector3d> allPositions = readPositionFile(result_path+"/positions_"+testName+".dat", NG*total_ns);
	vector<Vector4d> allRotations = readQuaternionFile(result_path+"/rotations_"+testName+".dat", NG*total_ns);

  {
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    // render window interacter
    // create renderer
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1.0, 1.0, 1.0); // white
    renderWindow->AddRenderer(renderer);

    std::mt19937 gen(123123);
  	std::uniform_real_distribution<double> rand_real(0., 1.0);

    vector<Vector3d> positions(allPositions.begin()+ns*NG, allPositions.begin()+(ns+1)*NG);
    vector<Vector4d> rotations(allRotations.begin()+ns*NG, allRotations.begin()+(ns+1)*NG);

		double range = std::numeric_limits<double>::min();
		int id = 0;
		for(int i = 0; i < positions.size(); ++i){
			Vector3d pos = positions[i];
			double grainRange = sqrt(pos(0)*pos(0) + pos(1)*pos(1));
			id = range < grainRange ? i : id;
			range = range < grainRange ? grainRange : range;
		}
		cout<<"furthest grain is: "<<id<<" range: "<<range<<endl;

		double height = std::numeric_limits<double>::min();
		int hid = 0;
		for(int i = 0; i < positions.size(); ++i){
			Vector3d pos = positions[i];
			double grainHeight = pos(2);
			hid = height < grainHeight ? i : hid;
			height = height < grainHeight ? grainHeight : height;
		}
		cout<<"highest grain is: "<<hid<<" height: "<<height<<endl;

		int deleteCount = 0;
		for(int nz = 0; nz < N; nz++){
			for(int ny = 0; ny < N; ny++){
				for(int nx = 0; nx < N; nx++){
					int disp_i = (N*N*nz + N*ny + nx) * ng;
					for (size_t g = 0; g < ng; g+=1) {
						int gid = g + disp_i;
						Vector3d pos = positions[gid];
						if(g % 500 == 0) cout<<g<<endl;
						stringstream fname;
				    fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/Polyhedrons/poly_"<<g+1<<".dat";
				    Polyhedron poly = readPolyFile(fname.str());
			      vtkSmartPointer<vtkPolyData> trianglePolyData = constructVtkPoly(poly, pos, rotations[gid]/rotations[gid].norm());
			      // Create mapper and actor
			      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			      mapper->SetInputData(trianglePolyData);
			      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			      actor->SetMapper(mapper);
						actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
						actor->GetProperty()->SetOpacity(1.);
						renderer->AddActor(actor);
			    }
				}
			}
		}
		cout<<" deleted "<<deleteCount<<" grains."<<endl;



		{
		  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
			planeSource->SetPoint1(1200.,0.,0.);
		  planeSource->SetPoint2(0.,1200.,0.);
		  planeSource->SetCenter(600.,600.,0.);
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
		  planeSource->SetPoint1(1200.,0.,0.);
		  planeSource->SetPoint2(0.,0.,1200.);
		  planeSource->SetCenter(600.,0.,600.);
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
		  planeSource->SetPoint1(0.,1200.,0.);
		  planeSource->SetPoint2(0.,0.,1200.);
		  planeSource->SetCenter(0.,600.,600.);
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
	  camera->SetPosition(100,100,150);
		//camera->SetPosition(0,0,200);
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
  	stringstream fname; fname << "/home/hasitha/Desktop/data/fabric/"+testName+"/images/1600steps_up2_"<<ns<<".png";

		//writer for saving later
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
  	writer->SetFileName(fname.str().c_str());
  	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
		writer->SetCompressionLevel(0);
  	writer->Write();

		vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	  renderWindowInteractor->SetRenderWindow(renderWindow);
	  renderWindowInteractor->Start();
  }

	cout << "program terminated" << endl;
	return EXIT_SUCCESS;
}
