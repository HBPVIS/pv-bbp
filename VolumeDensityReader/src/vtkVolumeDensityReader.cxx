
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkPointData.h"

// Header of the Reader
#include "vtkVolumeDensityReader.h"

#include <BBP/common.h>
#include <BBP/microcircuit.h>
#include <BBP/Filters/Microcircuit/Segment_Voxel_Mapper.h>
#include <BBP/Filters/Microcircuit/Segment_Voxel_Density.h>

using namespace bbp;

vtkCxxRevisionMacro(vtkVolumeDensityReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkVolumeDensityReader);

// ----------------------------------------------------------------------------
 
vtkVolumeDensityReader::vtkVolumeDensityReader()
{
    this->FileName = NULL;
   
    this->cx   = 419.0;   // C080301B1 : 419.134, 505.664, 1071.74
    this->cy   = 505.0;
    this->cz   = 1071.0;
    this->dimx = 50.0;
    this->dimy = 120.0;
    this->dimz = 50.0;
    this->res  = 1.0;

    // Set the whole extent
    this->WholeExtent[0] = 0;
    this->WholeExtent[1] = dimx / res;
    this->WholeExtent[2] = 0;
    this->WholeExtent[3] = dimy / res;
    this->WholeExtent[4] = 0;
    this->WholeExtent[5] = dimz / res;

    //this->SectionType = NULL;
    this->Target = NULL;

    this->SectionType = (char*)"all"; 

    this->SetTarget("viz1");

    this->SetNumberOfInputPorts(0);
}

// ----------------------------------------------------------------------------
 
vtkVolumeDensityReader::~vtkVolumeDensityReader()
{
}

//// ----------------------------------------------------------------------------
//
//void vtkVolumeDensityReader::SetWholeExtent(
//                const long unsigned int& min_x, int max_x,
//                const long unsigned int& min_y, int max_y,
//                const long unsigned int& min_z, int max_z)
//{
//    this->WholeExtent[0] = min_x;  this->WholeExtent[1] = max_x;
//    this->WholeExtent[2] = min_y;  this->WholeExtent[3] = max_y;
//    this->WholeExtent[4] = min_z;  this->WholeExtent[5] = max_z;
//}
//
//// ----------------------------------------------------------------------------
//
//void vtkVolumeDensityReader::SetSegmentVoxelMapper(float cx,   float cy,   float cz,
//                                                   float dimx, float dimy, float dimz,
//                                                   float res)
//{
//    this->cx   = cx;
//    this->cy   = cy;
//    this->cz   = cz;
//    this->dimx = dimx;
//    this->dimy = dimy;
//    this->dimz = dimz;
//    this->res  = res;
//    // Set the whole extent
//    this->WholeExtent[0] = 0;
//    this->WholeExtent[1] = dimx / res;
//    this->WholeExtent[2] = 0;
//    this->WholeExtent[3] = dimy / res;
//    this->WholeExtent[4] = 0;
//    this->WholeExtent[5] = dimz / res;
//
//}


// ----------------------------------------------------------------------------

int vtkVolumeDensityReader::RequestInformation(
   vtkInformation *vtkNotUsed(request),
   vtkInformationVector **vtkNotUsed(inputVector),
   vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
  
    int tmpExt[6], i;
    for (i = 0; i < 3; i++)
      {
      tmpExt[2*i] = this->WholeExtent[2*i];  
      tmpExt[2*i+1] = this->WholeExtent[2*i+1];
      }
  
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
                 tmpExt,6);
  
    vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_UNSIGNED_CHAR, 1);
  
    return 1;
}

// ----------------------------------------------------------------------------

void vtkVolumeDensityReader::ExecuteData(vtkDataObject *output)
{
    
    Timer timer;
    Pipeline pipeline;

    // open the Experiment file
    if (!this->FileName)
    {
        vtkErrorMacro(<< "A FileName must be specified.");
        return ;
    }
   
    const URI blue_config_file = this->FileName;
    Experiment experiment;
    experiment.open(blue_config_file);
    Microcircuit & microcircuit = experiment.microcircuit();

    std::cout << "Target = " << this->Target << std::endl;
    Cell_Target cell_target = experiment.cell_target(this->Target);
    timer.start("Loading microcircuit");
    microcircuit.load(cell_target, NEURONS | MORPHOLOGIES);
    timer.print();

    // define volume dimensions and resolution
    Box <Micron> bounding_box (this->cx, this->cy, this->cz, 
                               this->dimx, this->dimy, this->dimz);

    // define the section_type
    const char* section_type_string = this->SectionType;
    int sec_type = 0;
    if (!strcmp(section_type_string,"all") || !strcmp(section_type_string,"soma"))
        sec_type |= SOMA;
    if (!strcmp(section_type_string,"all") || !strcmp(section_type_string,"axons"))
        sec_type |= AXON;
    if (!strcmp(section_type_string,"all") || !strcmp(section_type_string,"dendrites"))
        sec_type |= DENDRITE | APICAL_DENDRITE;

    // configure and run processing pipeline
    timer.start("Setting up pipeline");
    bool clean_cut = true;
    Segment_Voxel_Mapper volume_mapping(this->res, bounding_box, clean_cut, Section_Type(sec_type));
    Segment_Voxel_Density <unsigned char> density_volume;
    volume_mapping.input_ptr()->reset(experiment.microcircuit_ptr());

    pipeline = volume_mapping + density_volume;// + nrrd_file_writer;
    timer.print();

    timer.start("Running pipeline");
    pipeline.mode(PROCESS_SINGLE_STEP_TRIGGER);
    pipeline.start();
    
    // Process the SDK-Pipeline 
    
    pipeline.step();
    timer.print();


    // Fill the vtkImageData with the Volume generated by the SDK Pipeline 
    timer.start("Filling VTK data structure");

    bbp::Voxels<unsigned char, Micron, Millisecond, TENSOR_LINEAR_VECTOR> volume = density_volume.output();
    
    int *whlExt;
    whlExt = this->GetWholeExtent();

    // Build the VTK image
    vtkImageData *data;
    data = this->AllocateOutputData(output);
    if (data->GetScalarType() != VTK_UNSIGNED_CHAR)
    {
      vtkErrorMacro("Execute: This source only outputs unsigned char");
      return;
    }
    if (data->GetNumberOfPoints() <= 0)
      return;
    
    data->GetPointData()->GetScalars()->SetName("MyData");

    for (int x=whlExt[0]; x<whlExt[1]; x++)
    {
        for (int y=whlExt[2]; y<whlExt[3]; y++)
        {
            for (int z=whlExt[4]; z<whlExt[5]; z++)
            {
                unsigned char* pixel = static_cast<unsigned char*>(data->GetScalarPointer(x,y,z));
                pixel[0] = volume(x,y,z);
            }
        }
    }
    timer.print();

    timer.start("Stopping pipeline");
    pipeline.stop();
    timer.print();

}
 
//----------------------------------------------------------------------------

void vtkVolumeDensityReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
 
  os << indent << "File Name: " 
      << (this->FileName ? this->FileName : "(none)") << "\n";  
}
