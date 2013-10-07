#include "vtkImageAlgorithmFilter.h"
 
#include "vtkImageData.h"

#include "vtkPointData.h" 
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
 
vtkCxxRevisionMacro(vtkImageAlgorithmFilter, "$Revision: 1.70 $");
vtkStandardNewMacro(vtkImageAlgorithmFilter);
 
vtkImageAlgorithmFilter::vtkImageAlgorithmFilter()
{
    this->WholeExtent[0] = -10.0;  this->WholeExtent[1] = 10.0;
    this->WholeExtent[2] = -10.0;  this->WholeExtent[3] = 10.0;
    this->WholeExtent[4] = -10.0;  this->WholeExtent[5] = 10.0;

    this->SetNumberOfInputPorts(0);

}
 
vtkImageAlgorithmFilter::~vtkImageAlgorithmFilter()
{
 
}

// ----------------------------------------------------------------------------

int vtkImageAlgorithmFilter::RequestInformation(
   vtkInformation *vtkNotUsed(request),
   vtkInformationVector **vtkNotUsed(inputVector),
   vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int tmpExt[6], i;
  for (i = 0; i < 3; i++)
    {
    tmpExt[2*i] = this->WholeExtent[2*i];  
    tmpExt[2*i+1] = this->WholeExtent[2*i+1];
    }

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               tmpExt,6);

  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 1);
  return 1;
}


// ----------------------------------------------------------------------------
 
void vtkImageAlgorithmFilter::ExecuteData(vtkDataObject *output)
{
    vtkImageData *data;
    int *whlExt;
    
    data = this->AllocateOutputData(output);
    if (data->GetScalarType() != VTK_FLOAT)
    {
      vtkErrorMacro("Execute: This source only outputs floats");
      return;
    }
    if (data->GetNumberOfPoints() <= 0)
    {
      return;
    }
    
    whlExt = this->GetWholeExtent();
    data->GetPointData()->GetScalars()->SetName("MyData");
    
  for (int x=whlExt[0]; x<whlExt[1]; x++)
    {
    for (int y=whlExt[2]; y<whlExt[3]; y++)
      {
        for (int z=whlExt[4]; z<whlExt[5]; z++)
        {
            float* pixel = static_cast<float*>(data->GetScalarPointer(x,y,z));
          pixel[0] = x+y+z; 
        }
    }
    }
}

