#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
 
int main(int argc, char *argv[])
{
  //Create three points. We will join (Origin and P0) with a red line and (Origin and P1) with a green line
  double origin[3] = {0.0, 0.0, 0.0};
  double p0[3] = {1.0, 0.0, 0.0};
  double p1[3] = {0.0, 1.0, 0.0};
  double p2[3] = {0.0, 1.0, 2.0};
  double p3[3] = {1.0, 2.0, 3.0};
  double p4[3] = {0.0, 0.0, 0.0};
  double p5[3] = {1.0, 0.0, 0.0};
  double p6[3] = {0.0, 2.0, 0.0};
  double p7[3] = {0.0, 0.0, 3.0};
 
  //create a vtkPoints object and store the points in it
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(origin);
  points->InsertNextPoint(p0);
  points->InsertNextPoint(p1);
  points->InsertNextPoint(p2);
  points->InsertNextPoint(p3);
  points->InsertNextPoint(p4);
  points->InsertNextPoint(p5);
  points->InsertNextPoint(p6);
  points->InsertNextPoint(p7);
 
  vtkSmartPointer<vtkPolyLine> polyLine1 = 
      vtkSmartPointer<vtkPolyLine>::New();
  polyLine1->GetPointIds()->SetNumberOfIds(5);
  for(unsigned int i = 0; i < 5; i++)
    {
    polyLine1->GetPointIds()->SetId(i,i);
    }
 
  vtkSmartPointer<vtkPolyLine> polyLine2 = 
      vtkSmartPointer<vtkPolyLine>::New();
  polyLine2->GetPointIds()->SetNumberOfIds(4);
  for(unsigned int i = 5; i < 9; i++)
    {
    polyLine2->GetPointIds()->SetId(i,i);
    }

  //Create a cell array to store the lines in and add the lines to it
  vtkSmartPointer<vtkCellArray> cells = 
      vtkSmartPointer<vtkCellArray>::New();
  cells->InsertNextCell(polyLine1);
  cells->InsertNextCell(polyLine2);
 
  //Create a polydata to store everything in
  vtkSmartPointer<vtkPolyData> polyData = 
      vtkSmartPointer<vtkPolyData>::New();
 
  //add the points to the dataset
  polyData->SetPoints(points);
 
  //add the lines to the dataset
  polyData->SetLines(cells);
 
  //setup actor and mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper = 
      vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(polyData);
 
  vtkSmartPointer<vtkActor> actor = 
      vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
 
  //setup render window, renderer, and interactor
  vtkSmartPointer<vtkRenderer> renderer = 
      vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = 
      vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderer->AddActor(actor);
 
  renderWindow->Render();
  renderWindowInteractor->Start();
 
  return EXIT_SUCCESS;
}
