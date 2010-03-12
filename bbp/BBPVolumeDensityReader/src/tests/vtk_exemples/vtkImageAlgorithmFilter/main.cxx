#include <vtkSmartPointer.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkOutlineFilter.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSampleFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkVolumeTextureMapper3D.h>
#include <vtkVolume.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkImageCast.h>

 
#include "vtkImageAlgorithmFilter.h"
#include "vtkImageData.h"
 
void PrintImage(vtkImageData* image);
 
int main (int argc, char *argv[])
{
  vtkSmartPointer<vtkImageAlgorithmFilter> filter = vtkSmartPointer<vtkImageAlgorithmFilter>::New();
  //filter->SetInput(input);
  filter->Update();
 
  vtkImageData* output = filter->GetOutput();
 
  vtkstd::cout << "Output image: " << vtkstd::endl;
  PrintImage(output);

        // bounding box  
    vtkSmartPointer<vtkOutlineFilter> outline = vtkSmartPointer<vtkOutlineFilter>::New();  
    outline->SetInput( output );  
    vtkSmartPointer<vtkPolyDataMapper> outlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();  
    outlineMapper->SetInput( outline->GetOutput() );
    vtkSmartPointer<vtkActor> outlineActor = vtkSmartPointer<vtkActor>::New();
    outlineActor->SetMapper( outlineMapper );

   // volume
    vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction
        = vtkSmartPointer<vtkPiecewiseFunction>::New();
       opacityTransferFunction->AddPoint(0, 0.0);
       opacityTransferFunction->AddPoint(5, 1.0);
   
    vtkSmartPointer<vtkColorTransferFunction> colorTansfertFunction
        = vtkSmartPointer<vtkColorTransferFunction>::New();
       colorTansfertFunction->AddRGBPoint(0, 1, 1, 1);
       colorTansfertFunction->AddRGBPoint(126, .665, .84, .88);
       colorTansfertFunction->AddRGBPoint(127, 0, 1, 0);
       colorTansfertFunction->AddRGBPoint(199, 1, 0, 0);
   
    vtkSmartPointer<vtkVolumeProperty> volumeProperty
        = vtkSmartPointer<vtkVolumeProperty>::New();
       volumeProperty->SetScalarOpacity(opacityTransferFunction);
       volumeProperty->SetColor(colorTansfertFunction);
       //volumeProperty->IndependentComponentsOff();
       volumeProperty->SetInterpolationTypeToLinear();
       volumeProperty->ShadeOn();

    vtkSmartPointer<vtkVolumeTextureMapper3D> volumeMapper
        = vtkSmartPointer<vtkVolumeTextureMapper3D>::New();
       volumeMapper->SetInput(output);

    vtkSmartPointer<vtkVolume> volume
        = vtkVolume::New();
       volume->SetMapper(volumeMapper);
       volume->SetProperty(volumeProperty);



      // a renderer window
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // an interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // add the actors to the scene
  renderer->AddActor(outlineActor);
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
 
//void PrintImage(vtkImageData* image)
//{
//  int* dims = image->GetDimensions();
// 
//  for (int y=0; y<dims[1]; y++)
//    {
//    for (int x=0; x<dims[0]; x++)
//      {
//      double v = image->GetScalarComponentAsDouble(x,y,0,0);
//      vtkstd::cout << v << " ";
//    }
//    vtkstd::cout << vtkstd::endl;
//  }
// 
//}

void PrintImage(vtkImageData* image)
{
  int* dims = image->GetDimensions();
 
  for (int y=-10; y<10; y++)
    {
    for (int x=-10; x<10; x++)
      {
        for (int z=-10; z<10; z++)
        {
            double v = image->GetScalarComponentAsDouble(x,y,z,0);
            vtkstd::cout << v << " ";
        }
    }
    vtkstd::cout << vtkstd::endl;
  }
 
}


