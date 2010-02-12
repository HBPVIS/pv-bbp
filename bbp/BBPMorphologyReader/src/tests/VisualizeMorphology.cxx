#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkTubeFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include "vtkMorphologyReader.h"

int main (int argc, char *argv[])
{

  //parse command line arguments
  if(argc != 2)
  {
  vtkstd::cout << "Required arguments: Filename" << vtkstd::endl;
  exit(-1);
  }
    
  vtkstd::string filename = argv[1];

  vtkSmartPointer<vtkMorphologyReader> reader = vtkSmartPointer<vtkMorphologyReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
 
  // Get the polyData
  vtkPolyData* polyData = reader->GetOutput();
  vtkSmartPointer<vtkTubeFilter> morph = vtkSmartPointer<vtkTubeFilter>::New();
  morph->SetInput(polyData);
  morph->SetNumberOfSides(8);
  morph->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
 
  //create a mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(morph->GetOutputPort());
  mapper->ScalarVisibilityOn();
  mapper->SetScalarModeToUsePointFieldData();
  mapper->SelectColorArray("Colors");

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
  renderer->SetBackground(.2, .3, .4); 
 
  // render an image (lights and cameras are created automatically)
  renderWindow->Render();

  // change the interactor style
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = 
    vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(style); 

  // begin mouse interaction
  renderWindowInteractor->Start();
 
  return 0;
}
