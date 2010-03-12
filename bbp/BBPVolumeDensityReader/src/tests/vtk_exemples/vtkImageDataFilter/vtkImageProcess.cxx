
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"

// Header of the Reader
#include "vtkImageProcess.h"

vtkCxxRevisionMacro(vtkImageProcess, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkImageProcess);
 
vtkImageProcess::vtkImageProcess()
{
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

vtkImageProcess::~vtkImageProcess()
{
}


int vtkImageProcess::FillOutputPortInformation( int port, vtkInformation* info )
{
  if ( port == 0 )
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData" );
 
    return 1;
  }
 
  return 0;
}

int vtkImageProcess::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
 
    // get the ouptut
    vtkImageData *output = vtkImageData::SafeDownCast(
            outInfo->Get(vtkDataObject::DATA_OBJECT()));
    
    vtkSmartPointer < vtkImageData > imageData = vtkSmartPointer<vtkImageData>::New();
    //specify the size of the image data
    imageData->SetDimensions(50,100,50);
    imageData->SetNumberOfScalarComponents(1);
    imageData->SetScalarTypeToUnsignedChar();
  
    int* dims = imageData->GetDimensions();
  
    cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << endl;
    cout << "Number of points: " << imageData->GetNumberOfPoints() << endl;
    cout << "Number of cells: " << imageData->GetNumberOfCells() << endl;
  
    for (int z = 0; z < dims[2]; z++)
      {
        for (int y = 0; y < dims[1]; y++)
        {
        for (int x = 0; x < dims[0]; x++)
          {
          unsigned char* pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(x,y,z));
          pixel[0] = 'd';
          }
        }
      }

    //output->ShallowCopy(imageData);
    output->DeepCopy(imageData);


    return 1;
}
 
 
//----------------------------------------------------------------------------
void vtkImageProcess::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
