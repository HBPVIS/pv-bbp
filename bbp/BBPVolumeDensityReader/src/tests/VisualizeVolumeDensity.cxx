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

#include "vtkVolumeDensityReader.h"

int main (int argc, char *argv[])
{
    //parse command line arguments
    if(argc != 2)
    {
    vtkstd::cout << "Required arguments: Filename" << vtkstd::endl;
    exit(-1);
    }

    vtkstd::string filename = argv[1];

    int* dims;

    vtkSmartPointer<vtkVolumeDensityReader> reader = vtkSmartPointer<vtkVolumeDensityReader>::New();
    reader->SetFileName(filename.c_str());
    //reader->SetSegmentVoxelMapper(275.0,600.0,275.0,500.0,1200.0,500.0,2.0);
    reader->SetSectionType("all");
    reader->SetTarget("Layer5");
    reader->SetTarget("MiniColum_0");
    reader->Update();

    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData = reader->GetOutput();
    dims = imageData->GetDimensions();

    imageData->SetScalarTypeToUnsignedChar();

    //cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << endl;
    //cout << "Number of points: " << imageData->GetNumberOfPoints() << endl;
    //cout << "Number of cells: " << imageData->GetNumberOfCells() << endl;

    // bounding box  
    vtkSmartPointer<vtkOutlineFilter> outline = vtkSmartPointer<vtkOutlineFilter>::New();
    outline->SetInput( imageData );
    vtkSmartPointer<vtkPolyDataMapper> outlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    outlineMapper->SetInput( outline->GetOutput() );
    vtkSmartPointer<vtkActor> outlineActor = vtkSmartPointer<vtkActor>::New();
    outlineActor->SetMapper( outlineMapper );

    // volume
    vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction
        = vtkSmartPointer<vtkPiecewiseFunction>::New();
       opacityTransferFunction->AddPoint(0, 0.0);
       opacityTransferFunction->AddPoint(255, 1.0);

    vtkSmartPointer<vtkColorTransferFunction> colorTansfertFunction
        = vtkSmartPointer<vtkColorTransferFunction>::New();
       colorTansfertFunction->AddRGBPoint(0, 1, 1, 1);
       colorTansfertFunction->AddRGBPoint(255, 1, 1, 1);

   vtkSmartPointer<vtkVolumeProperty> volumeProperty
        = vtkSmartPointer<vtkVolumeProperty>::New();
       volumeProperty->SetScalarOpacity(opacityTransferFunction);
       volumeProperty->SetColor(colorTansfertFunction);
       //volumeProperty->IndependentComponentsOff();
       volumeProperty->SetInterpolationTypeToLinear();
       volumeProperty->ShadeOn();

    vtkSmartPointer<vtkVolumeTextureMapper3D> volumeMapper
        = vtkSmartPointer<vtkVolumeTextureMapper3D>::New();
       volumeMapper->SetInput(imageData);

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
