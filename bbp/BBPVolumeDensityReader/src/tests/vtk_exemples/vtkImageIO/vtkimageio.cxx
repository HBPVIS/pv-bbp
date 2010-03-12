#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkImageData.h"
#include "vtkFloatArray.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkPointData.h"
#include "vtkImageWriter.h"
#include "vtkImageReader.h"
#include "vtkSphereSource.h"
#include "vtkStructuredPointsWriter.h"
#include "vtkStructuredPointsReader.h"
#include "vtkDataSet.h"
#include "vtkStructuredData.h"
#include "vtkImageToStructuredPoints.h"
#include "vtkStructuredPointsSource.h"
#include "vtkDataArray.h"

main ()
{

float sp, z, y, x, s, point[3], point2[3];
int i, j, k,  kOffset, jOffset, offset,h,t=0;
double valuerange[2];

vtkImageData *vol = vtkImageData::New();
vol->SetDimensions(26,26,26);
vol->SetOrigin(-0.5,-0.5,-0.5);
sp =1.0/25.0;
vol->SetSpacing(sp,sp,sp);
vol->AllocateScalars();
vol->SetScalarTypeToFloat();

vtkFloatArray *scalars =vtkFloatArray::New();
for(k=0;k<26;k++)
{
z = -0.5 + k*sp;
kOffset = k*26*26;

for(j=0; j<26; j++)
{
y = -0.5 +j*sp;
jOffset = j*26;


for(i=0;i<26;i++)
{
x = -0.5 + i*sp;
s = x*x + y*y + z*z - (0.4*0.4);
offset = i + jOffset + kOffset;
scalars->InsertTuple1(offset,s);
//  scalars->InsertValue(offset,s);
      //   if(t<100)  printf("scalars -> %d , %lf\n",offset,s);
       // t++;
}
}
}

vol->GetPointData()->SetScalars(scalars);
vol->GetPointData()->CopyAllOn();
vol->Update();
vol->UpdateData();

vtkImageWriter *writer = vtkImageWriter::New();
writer->SetInput(vol);
writer->SetFileName("voldataim.vti");
writer->SetFileDimensionality(3);
writer->Write();

vtkImageReader *reader  = vtkImageReader::New();
reader->SetFileName("voldataim.vti");
reader->SetDataScalarTypeToFloat (); //very important !
reader->SetFileDimensionality(3);
reader->SetDataExtent(0,25, 0, 25, 0, 25);
reader->SetDataSpacing(1.0/25.0, 1.0/25.0, 1.0/25.0);
reader->SetDataOrigin(-0.5,-0.5,-0.5);

reader->GetOutput()->GetScalarRange(valuerange);//new

reader->Update();
reader->UpdateWholeExtent();
reader->UpdateInformation();

vtkFloatArray *scalars2 =vtkFloatArray::New();
scalars2=(vtkFloatArray *)reader->GetOutput()->GetPointData()->GetScalars();
for(i=0;i<100;i++) printf("reader->GetOutput())->GetPointData()->GetScalars(%d) = \
%lf\n", i, scalars2->GetValue(i)); 

vtkContourFilter *contour = vtkContourFilter::New();
contour->SetInput(reader->GetOutput());
//   contour->SetInput(vol);
contour->SetValue(0,0);
contour->Update();

vtkPolyDataMapper *volMapper =  vtkPolyDataMapper::New();
volMapper->SetInput(contour->GetOutput());
volMapper->ScalarVisibilityOff();
vtkActor *volActor = vtkActor::New();
volActor->SetMapper(volMapper);
// volActor->GetProperty()->SetRepresentationToWireframe();

vtkRenderer *renderer = vtkRenderer::New();
vtkRenderWindow *renWin = vtkRenderWindow::New();
renWin->AddRenderer(renderer);
vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
iren->SetRenderWindow(renWin);


renderer->AddActor(volActor);
renderer->SetBackground(1,1,1);
renWin->SetSize(400,400);

// interact with data

renWin->Render();
iren->Start();

// Clean up
renderer->Delete();
renWin->Delete();
iren->Delete();
vol->Delete();
scalars->Delete();
contour->Delete();
volMapper->Delete();
volActor->Delete();
}

