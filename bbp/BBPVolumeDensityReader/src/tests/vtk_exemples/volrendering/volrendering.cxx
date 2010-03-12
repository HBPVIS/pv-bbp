#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkOutlineFilter.h>
#include <vtkContourFilter.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkQuadric.h>
#include <vtkSampleFunction.h>
#include <vtkTexture.h>
#include <vtkPlaneSource.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPiecewiseFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkVolumeTextureMapper3D.h>
#include <vtkVolume.h>


int main()
{


////Create image data
//vtkImageData* image = vtkImageData::New();
//   image->SetDimensions(sizeX, sizeY, sizeZ);
//   image->SetOrigin(0,0,0);
//   image->SetSpacing(0.1,0.1,0.1);
//   image->SetNumberOfScalarComponents(4);
//   image->SetScalarTypeToUnsignedChar();
//   image->GetPointData()->SetScalars(createData());

  // create the quadric function definition
  vtkQuadric *quadric = vtkQuadric::New();
  quadric->SetCoefficients(.5,1,.2,0,.1,0,0,.2,0,0);

  // sample the quadric function
  vtkSampleFunction *sample = vtkSampleFunction::New();
  sample->SetSampleDimensions(50,50,50);
  sample->SetImplicitFunction(quadric);

  vtkImageData* imageData = sample->GetOutput();



vtkPiecewiseFunction* opacityTransferFunction = vtkPiecewiseFunction::New();
   opacityTransferFunction->AddPoint(0, 0.0);
   opacityTransferFunction->AddPoint(5, 1.0);


vtkVolumeProperty* volumeProperty = vtkVolumeProperty::New();
   volumeProperty->SetScalarOpacity(opacityTransferFunction);
   volumeProperty->IndependentComponentsOff();
      

vtkVolumeTextureMapper3D* volumeMapper = vtkVolumeTextureMapper3D::New();
   volumeMapper->SetInput(imageData);        

vtkVolume* volume = vtkVolume::New();
   volume->SetMapper(volumeMapper);
   volume->SetProperty(volumeProperty); 

  // a renderer and render window
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // an interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // add the actors to the scene
  renderer->AddVolume(volume);
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
