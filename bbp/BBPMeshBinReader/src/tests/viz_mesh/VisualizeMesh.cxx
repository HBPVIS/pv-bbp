#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>

#include "vtkMeshBinReader.h"

int main (int argc, char *argv[])
{

  //parse command line arguments
  if(argc != 2)
  {
  vtkstd::cout << "Required arguments: Filename" << vtkstd::endl;
  exit(-1);
  }
    
  vtkstd::string filename = argv[1];

  vtkSmartPointer<vtkMeshBinReader> reader = vtkSmartPointer<vtkMeshBinReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkPolyData* polydata = reader->GetOutput();

 
  //create a mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(polydata);
 
  // create an actor
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
 
  // a renderer and render window
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
 
  // an interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
 
  // add the actors to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(1,0.8,0.8); // Background color white
 
  // render an image (lights and cameras are created automatically)
  renderWindow->Render();
 
  // begin mouse interaction
  renderWindowInteractor->Start();
 
  return 0;
}
