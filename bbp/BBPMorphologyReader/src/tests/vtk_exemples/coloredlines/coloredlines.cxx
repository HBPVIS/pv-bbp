#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
 
int main(int, char *[])
{
  //Create three points. Join (Origin and P0) with a red line and
  //(Origin and P1) with a green line
  double origin[3] = {0.0, 0.0, 0.0};
  double p0[3] = {1.0, 0.0, 0.0};
  double p1[3] = {0.0, 1.0, 0.0};
 
  //Create a vtkPoints object and store the points in it
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  pts->InsertNextPoint(origin);
  pts->InsertNextPoint(p0);
  pts->InsertNextPoint(p1);
 
  //Setup two colors - one for each line
  unsigned char red[3] = {255, 0, 0};
  unsigned char green[3] = {0, 255, 0};
 
  //Setup the colors array
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");
 
  //Add the colors we created to the colors array
  colors->InsertNextTupleValue(red);
  colors->InsertNextTupleValue(green);
 
  //Create the first line (between Origin and P0)
  vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
  line0->GetPointIds()->SetId(0,0); //the second 0 is the index of the Origin in the vtkPoints
  line0->GetPointIds()->SetId(1,1); //the second 1 is the index of P0 in the vtkPoints
 
  //Create the second line (between Origin and P1)
  vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
  line1->GetPointIds()->SetId(0,0); //the second 0 is the index of the Origin in the vtkPoints
  line1->GetPointIds()->SetId(1,2); //2 is the index of P1 in the vtkPoints

      //Create the first line (between Origin and P0)
  vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
  line2->GetPointIds()->SetId(0,0); //the second 0 is the index of the Origin in the vtkPoints
  line2->GetPointIds()->SetId(1,1); //the second 1 is the index of P0 in the vtkPoints

  //Create the second line (between Origin and P1)
  vtkSmartPointer<vtkLine> line3 = vtkSmartPointer<vtkLine>::New();
  line3->GetPointIds()->SetId(0,0); //the second 0 is the index of the Origin in the vtkPoints
  line3->GetPointIds()->SetId(1,2); //2 is the index of P1 in the vtkPoints

  //Create a cell array to store the lines in and add the lines to it
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  lines->InsertNextCell(line0);
  lines->InsertNextCell(line1);
  //lines->InsertNextCell(line2);
  //lines->InsertNextCell(line3);
 
  //Create a polydata to store everything in
  vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
 
  //Add the points to the dataset
  linesPolyData->SetPoints(pts);
 
  //Add the lines to the dataset
  linesPolyData->SetLines(lines);
 
  //Color the lines - associate the first component (red) of the
  //colors array with the first component of the cell array (line 0)
  //and the second component (green) of the colors array with the
  //second component of the cell array (line 1)
  linesPolyData->GetCellData()->SetScalars(colors);
 
  //Setup actor and mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(linesPolyData);
 
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
 
  //Setup render window, renderer, and interactor
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderer->AddActor(actor);
 
  renderWindow->Render();
  renderWindowInteractor->Start();
 
  return EXIT_SUCCESS;
}
