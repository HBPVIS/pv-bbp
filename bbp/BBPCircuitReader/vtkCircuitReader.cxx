// Visual studio debug settings for paths specifi to this module
//
// PATH=D:\build\paraview\bin\Debug;C:\Program Files\hdf5-vfd-1.8.9\bin;%PATH%
// PV_PLUGIN_PATH=D:\build\buildyard\ParaBBP\bin\Debug
//

#include "vtkCellArray.h" 
#include "vtkTriangle.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyLine.h>
#include <vtkLine.h>
#include <vtkCellData.h>

//
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDataArraySelection.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"
//
#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkLongLongArray.h"
#include "vtkUnsignedLongLongArray.h"
#include "vtkIntArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkOutlineSource.h"
#include "vtkAppendPolyData.h"
#include "vtkBoundingBox.h"
//
#include <vtksys/SystemTools.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vtkstd/vector>
//
#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCharArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"
#include "vtkExtentTranslator.h"
#include "vtkBoundingBox.h"
#include "vtkIdListCollection.h"
//
#include "vtkDummyController.h"
//
#include <functional>
#include <algorithm>
#include <numeric>
#include <iterator>

// BBP-SDK Morphology Reader
#include "BBP/IO/File/Parsers/Morphology_HDF5_File_Parser.h"

// Header of the Reader
#include "vtkCircuitReader.h"

// BBP-SDK

#include <BBP/common.h>
#include "BBP/Microcircuit/Experiment.h"
#include "BBP/Microcircuit/Targets/Targets.h"
#include "BBP/Microcircuit/Targets/Cell_Target.h"
#include <BBP/Microcircuit/Containers/Neurons.h>

// Voxelization
#include "BBP/Voxelization/voxelization.h"
#include "BBP/VtkDebugging/visualization.h"


using namespace bbp;
//----------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkCircuitReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkCircuitReader);
//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkCircuitReader, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkCircuitReader::vtkCircuitReader()
{
  this->FileName = NULL;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  //
  this->NumberOfTimeSteps               = 0;
  this->TimeStep                        = 0;
  this->ActualTimeStep                  = 0;
  this->TimeStepTolerance               = 1E-6;
  this->FileName                        = NULL;
  this->UpdatePiece                     = 0;
  this->UpdateNumPieces                 = 0;
  this->IntegerTimeStepValues           = 0;
  this->PointDataArraySelection         = vtkDataArraySelection::New();
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
}
//----------------------------------------------------------------------------
vtkCircuitReader::~vtkCircuitReader()
{
  delete []this->FileName;
  this->SetController(NULL);
}
//----------------------------------------------------------------------------
int vtkCircuitReader::FillOutputPortInformation( int port, vtkInformation* info )
{
  if ( port == 0 )
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );

    return 1;
  }

  return 0;
} 
//----------------------------------------------------------------------------
int vtkCircuitReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{

  // if there is no blue config supplied yet, exit quietly
  if (!this->FileName)
  {
    return 1;
  }

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  bool NeedToReadInformation = (FileModifiedTime>FileOpenedTime);

  if (NeedToReadInformation)
  {

    // ----------------------------------------------------------------------
    // Set parameters
    // ----------------------------------------------------------------------

    std::string                 blueconfig           = this->FileName;
    std::string                 target_name          = "AllCompartments";

    // -------------------------------------------------------------------   
    // Create BBP-SDK Experiment and Microcircuit to access to the neurons.
    // -------------------------------------------------------------------   
    experiment.open(blueconfig);
    bbp::Microcircuit& microcircuit = experiment.microcircuit();
    bbp::Target target = experiment.user_targets().get_target(target_name);
    microcircuit.load(target, bbp::NEURONS | bbp::MORPHOLOGIES | bbp::MESHES);
    bbp::Neurons & neurons = microcircuit.neurons(); 


    this->NumberOfTimeSteps = 1;
    this->TimeStepValues.assign(this->NumberOfTimeSteps, 0.0);
    int validTimes = 0;
    for (int i=0; i<this->NumberOfTimeSteps; ++i)
    {
    }

    if (this->NumberOfTimeSteps==0)
    {
      vtkErrorMacro(<<"No time steps in data");
      return 0;
    }

    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
      &this->TimeStepValues[0],
      static_cast<int>(this->TimeStepValues.size()));
    double timeRange[2];
    timeRange[0] = this->TimeStepValues.front();
    timeRange[1] = this->TimeStepValues.back();
    if (this->TimeStepValues.size()>1)
    {
      this->TimeStepTolerance = 0.01*(this->TimeStepValues[1]-this->TimeStepValues[0]);
    }
    else
    {
      this->TimeStepTolerance = 1E-3;
    }
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkCircuitReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{

  // open the HDF5 file
  if (!this->FileName)
  {
    vtkErrorMacro(<< "A FileName must be specified.");
    return 0;
  }

  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the ouptut
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  bbp::Microcircuit& microcircuit = experiment.microcircuit();
  bbp::Neurons & neurons = microcircuit.neurons(); 

  // Set VTK points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkIntArray> neuronId = vtkSmartPointer<vtkIntArray>::New();
  neuronId->SetName("NeuronId");
  // Set VTK triangles
  vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

  std::cout << "Neuron count : " << neurons.size() << std::endl;

  // -------------------------------------------------------------------   
  // iterate over neurons
  // -------------------------------------------------------------------   
  int Ncount = 0;
  for (bbp::Neurons::iterator i = neurons.begin(); i != neurons.end(); ++i,++Ncount)
  {
    if (Ncount>25) break;
    //
    const bbp::Morphology* sdk_morph    = &i->morphology();
    const bbp::Sections    sdk_sections = i->neurites();
    const bbp::Mesh*       sdk_mesh     = &i->morphology().mesh();
    Degree_Angle           y_rotation   = i->orientation().rotation * (-1.0);
    const Vector3 position( i->position().x(), i->position().y(), i->position().z( ));

    bbp::Vertex_Index vertexCount = sdk_mesh->vertex_count();
    bbp::Triangle_Index faceCount = sdk_mesh->triangle_count();
    bbp::Triangle_Index stripLength = sdk_mesh->triangle_strip_length();
    //   Vector_3D_Micron_Array vertices;
    //   Section_ID_Array vertex_sections;
    //   Float_Array vertex_relative_distances;
    //   Vertex_Index_Array faces;
    //   Vertex_Index_Array strips;

    const Vector_3D<bbp::Micron>* vertices2 = sdk_mesh->vertices().pointer(); 

    std::cout << "Neuron : " << std::distance(neurons.begin(),i) << std::endl;
    std::cout << "vertex count : " << vertexCount << std::endl;
    vtkIdType oldN =  points->GetNumberOfPoints();
    vtkIdType newN =  oldN + vertexCount;
    //
    points->GetData()->Resize(newN);
    points->SetNumberOfPoints(newN);
    //
    neuronId->Resize(newN);
    neuronId->SetNumberOfTuples(newN);
    //
    vtkstd::cout  << faceCount << std::endl;
    //
    vtkIdType Ncell =  triangles->GetNumberOfCells();
    vtkIdType *cells = triangles->WritePointer(Ncell+faceCount, 4*(Ncell+faceCount));
    //
    vtkIdType ip = oldN;
    for ( bbp::Vertex_Index v = 0 ; v < vertexCount ; ++v )
    {
      points->SetPoint(ip, vertices2[v].x(), vertices2[v].y(), vertices2[v].z());
      neuronId->SetValue(ip,Ncount);
      ip++;
    }

    const bbp::Vertex_Index* faces2 = sdk_mesh->triangles().pointer();
    for ( bbp::Triangle_Index t = 0 ; t < faceCount ; ++t )
    {
      cells[Ncell*4 + 0] = 3;
      cells[Ncell*4 + 1] = oldN + faces2[t*3 + 0];
      cells[Ncell*4 + 2] = oldN + faces2[t*3 + 1];
      cells[Ncell*4 + 3] = oldN + faces2[t*3 + 2];
      Ncell++;
    }  
  }
  output->SetPolys(triangles);

/*
  vtkIdType Np =  points->GetNumberOfPoints();
  {
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType *cells = vertices->WritePointer(Np, 2*Np);
    for (vtkIdType i=0; i<Np; ++i)
    {
      cells[2*i] = 1;
      cells[2*i+1] = i;
    }
    output->SetVerts(vertices);
  }
*/

  // Set vtkPolyData
  //  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  output->SetPoints(points);
  output->GetPointData()->AddArray(neuronId);
  //  polydata->SetPolys(triangles);

  //output = polydata; //doesn't work
  //  output->ShallowCopy(polydata);

  return 1;
}


//----------------------------------------------------------------------------
void vtkCircuitReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
    << (this->FileName ? this->FileName : "(none)") << "\n";  
}
//----------------------------------------------------------------------------
void vtkCircuitReader::SetFileName(char *filename)
{
  if (this->FileName == NULL && filename == NULL)
  {
    return;
  }
  if (this->FileName && filename && (!strcmp(this->FileName,filename)))
  {
    return;
  }
  delete [] this->FileName;
  this->FileName = NULL;

  if (filename)
  {
    this->FileName = vtksys::SystemTools::DuplicateString(filename);
    this->SetFileModified();
  }
  this->Modified();
}
//----------------------------------------------------------------------------
void vtkCircuitReader::SetFileModified()
{
  this->FileModifiedTime.Modified();
  this->Modified();
}

//----------------------------------------------------------------------------
const char* vtkCircuitReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkCircuitReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkCircuitReader::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name))
  {
    if (status)
    {
      this->PointDataArraySelection->EnableArray(name);
    }
    else
    {
      this->PointDataArraySelection->DisableArray(name);
    }
    this->Modified();
  }
}
//----------------------------------------------------------------------------
void vtkCircuitReader::Enable(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}
//----------------------------------------------------------------------------
void vtkCircuitReader::Disable(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}
//----------------------------------------------------------------------------
void vtkCircuitReader::EnableAll()
{
  this->PointDataArraySelection->EnableAllArrays();
}
//----------------------------------------------------------------------------
void vtkCircuitReader::DisableAll()
{
  this->PointDataArraySelection->DisableAllArrays();
}
//----------------------------------------------------------------------------
int vtkCircuitReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
