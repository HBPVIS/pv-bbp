// Visual studio debug settings for paths specific to this module
//
// PATH=D:\build\paraview-3.98\bin\Debug;C:\Program Files\hdf5-vfd-1.8.9\bin;%PATH%
// PV_PLUGIN_PATH=D:\build\buildyard\ParaBBP\bin\Debug
// _NO_DEBUG_HEAP=1
// working directory : D:\build\paraview-3.98\bin\Debug
//

#include "vtkObjectFactory.h"
#include "vtkCellArray.h" 
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
//
#include "vtkDataArraySelection.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
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
#include "vtkVariantArray.h"
#include "vtkStringArray.h"
#include "vtkCellArray.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkExtentTranslator.h"
//
#include "vtkTransform.h"
#include "vtkPolyDataNormals.h"
#include "vtkDistributedDataFilter.h"
//
#include "vtkDummyController.h"
//
#include "vtkPKdTree.h"
#include "vtkBoundsExtentTranslator.h"
//
#include <vtksys/SystemTools.hxx>
//
#include <vector>
#include <deque>
#include <functional>
#include <algorithm>
#include <numeric>
#include <iterator>
//
// BBP-SDK Morphology Reader
#include "BBP/IO/File/Parsers/Morphology_HDF5_File_Parser.h"
#include "BBP/Microcircuit/Readers/Microcircuit_Reader.h"
#include "BBP/Microcircuit/Readers/Mesh_Reader.h"
#include "BBP/Microcircuit/Readers/compartmentReportReader.h"
#include "BBP/Microcircuit/Mappings/Compartment_Report_Mapping.h"
#include "BBP/Microcircuit/Soma.h"
#include "BBP/Microcircuit/Mesh.h"
#include "BBP/Microcircuit/Datasets/Morphology_Dataset.h"

// Header of this Reader
#include "vtkCircuitReader.h"
//#include "SpikeData.h"

// BBP-SDK
// Voxelization
//#include "BBP/Voxelization/voxelization.h"
//#include "BBP/VtkDebugging/visualization.h"

#define MANUAL_MESH_LOAD

//----------------------------------------------------------------------------
#define BBP_ARRAY_NAME_NORMAL           "Normal"
#define BBP_ARRAY_NAME_NEURONGID        "Neuron Gid"
#define BBP_ARRAY_NAME_NEURONINDEX      "Neuron Index"
#define BBP_ARRAY_NAME_SECTION_ID       "Section Id"
#define BBP_ARRAY_NAME_SECTION_TYPE     "Section Type"
#define BBP_ARRAY_NAME_DENDRITE_RADIUS  "Dendrite Radius"
#define BBP_ARRAY_NAME_VOLTAGE          "Voltage"
#define BBP_ARRAY_NAME_RTNEURON_OPACITY "RTNeuron Opacity"

//----------------------------------------------------------------------------
#define USE_BBP_TRANSFORM
#define USE_VTK_TRANSFORM
//----------------------------------------------------------------------------
using namespace bbp;
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkCircuitReader);
//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkCircuitReader, Controller, vtkMultiProcessController);
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
  this->CurrentTime                     = std::numeric_limits<float>::min();
  this->TimeStepTolerance               = 1E-6;
  this->FileName                        = NULL;
  this->DefaultTarget                   = NULL;
  this->UpdatePiece                     = 0;
  this->UpdateNumPieces                 = 0;
  this->IntegerTimeStepValues           = 0;
  this->ExportNeuronMesh                = 1;
  this->ExportMorphologySkeleton        = 0;
  this->ParallelRedistribution          = 1;
  this->MaximumNumberOfNeurons          = 25;
  //
  this->PointDataArraySelection         = vtkDataArraySelection::New();
  this->TargetsSelection                = vtkDataArraySelection::New();
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
  this->SIL              = vtkSmartPointer<vtkMutableDirectedGraph>::New();
  this->SILUpdateStamp   = 0;
  //
  this->CachedNeuronMesh              = vtkSmartPointer<vtkPolyData>::New();
  this->CachedMorphologySkeleton      = vtkSmartPointer<vtkPolyData>::New();
  this->DistributedDataFilter         = NULL;
  this->BoundsTranslator              = vtkSmartPointer<vtkBoundsExtentTranslator>::New();
}
//----------------------------------------------------------------------------
vtkCircuitReader::~vtkCircuitReader()
{
  this->SIL                           = NULL;
  this->CachedNeuronMesh              = NULL;
  this->CachedMorphologySkeleton      = NULL;
  this->DistributedDataFilter         = NULL;
  //
  this->PointDataArraySelection->Delete();
  this->TargetsSelection->Delete();
  this->SetController(NULL);
  //
  delete []this->FileName;
  delete []this->DefaultTarget;
}
//----------------------------------------------------------------------------
int vtkCircuitReader::FillOutputPortInformation( int port, vtkInformation* info )
{
  if (port == 0) {
    if (this->Controller && this->Controller->GetNumberOfProcesses()>1) {
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
    }
    else {
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
    }
    return 1;
  }

  if (port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
    return 1;
  }

  return 0;
} 
/*
std::pair<double, double> GetClosest(std::vector<double> &sortedlist, const double& val) const
{
  std::vector<double>::const_iterator it = std::lower_bound(sortedlist.begin(), sortedlist.end(), val);
  if (it == sortedlist.end())        return std::make_pair(sortedlist.back(), sortedlist.back());
  else if (it == sortedlist.begin()) return std::make_pair(sortedlist.front(), sortedlist.front());
  else return std::make_pair(*(it - 1), *(it));
}
*/
class vtkCircuitReaderToleranceCheck: public std::binary_function<double, double, bool>
{
public:
  vtkCircuitReaderToleranceCheck(double tol) { this->tolerance = tol; }
  double tolerance;
  //
    result_type operator()(first_argument_type a, second_argument_type b) const
    {
      bool result = (fabs(a-b)<=(this->tolerance));
      return (result_type)result;
    }
};
//----------------------------------------------------------------------------
int vtkCircuitReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  int result = 1;
  // if there is no blue config supplied yet, exit quietly
  if (!this->FileName) {
    return 1;
  }
  vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
//  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
  //
  this->UpdatePiece = this->Controller->GetLocalProcessId(); // outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = this->Controller->GetNumberOfProcesses(); // outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  //
  outInfo0->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  bool NeedToReadInformation = (FileModifiedTime>FileOpenedTime);

  if (!vtksys::SystemTools::FileExists(this->FileName)) {
    vtkWarningMacro("File not found " << this->FileName);
    NeedToReadInformation = 0;
    result = 0;
  }

  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_NORMAL);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_NEURONGID);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_NEURONINDEX);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_SECTION_ID);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_SECTION_TYPE);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_DENDRITE_RADIUS);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_VOLTAGE);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_RTNEURON_OPACITY);

  if (NeedToReadInformation) {
    std::string blueconfig = this->FileName;

    // -------------------------------------------------------------------   
    // Create BBP-SDK Experiment and Microcircuit to access to the neurons.
    // -------------------------------------------------------------------   
    this->Experiment.open(blueconfig);
    this->Microcircuit = this->Experiment.microcircuit_ptr();

    // default Target?
    this->TargetName = this->DefaultTarget ? this->DefaultTarget : this->TargetsSelection->GetArrayName(0);
    //
    int N = this->TargetsSelection->GetNumberOfArrays();
    for (int i=0; i<N; i++) {
      const char *name = this->TargetsSelection->GetArrayName(i);
      if (this->TargetsSelection->ArrayIsEnabled(name)) {
        this->TargetName = name;
        break;
        //std::cout << "Target selected : " << name << std::endl;
      }
    }

    try {
      this->Target = this->Experiment.user_targets().get_target(this->TargetName);
    }
    catch (...) {
      this->Target = this->Experiment.targets().get_target(this->TargetName);
    }

    // Don't load meshes yet, we'll do that once we've decided which neurons this node will generate
    this->Microcircuit->load(this->Target, 0); // bbp::NEURONS);

// std::cout <<"Made it past load " << std::endl;

    // time steps of reports are in the report file
    if (this->UpdateNumPieces==1) {
//      if (this->OpenReportFile()) {
//        this->NumberOfTimeSteps = (this->stopTime-this->startTime)/this->timestep;
//      }
//      else {
        this->NumberOfTimeSteps = 0;
//      }
    }
    else {
      this->NumberOfTimeSteps = 0;
    }

    if (this->NumberOfTimeSteps==0) {
//      vtkErrorMacro(<<"No time steps report data, may cause crash later : TODO fix this");
    }
    else {
      this->TimeStepValues.assign(this->NumberOfTimeSteps, 0.0);
      for (int i=0; i<this->NumberOfTimeSteps; ++i) {
        this->TimeStepValues[i] = this->startTime + (i*this->timestep);
      }
    }

    if (this->TimeStepValues.size()>0) {
      outInfo0->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
        &this->TimeStepValues[0],
        static_cast<int>(this->TimeStepValues.size()));
  //    outInfo1->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
  //      &this->TimeStepValues[0],
  //      static_cast<int>(this->TimeStepValues.size()));
      double timeRange[2] = { this->TimeStepValues.front(), this->TimeStepValues.back() };
      if (this->TimeStepValues.size()>1)
      {
        this->TimeStepTolerance = 0.01*(this->TimeStepValues[1]-this->TimeStepValues[0]);
      }
      else
      {
        this->TimeStepTolerance = 1E-3;
      }
      outInfo0->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    }
    //    outInfo1->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

    this->FileOpenedTime.Modified();
    result = 1;
  }

  if (this->Controller->GetNumberOfProcesses()>1 && this->ParallelRedistribution) {
    outInfo0->Set(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR(), this->BoundsTranslator);
  }

  //
//  this->BuildSIL();
//  outInfo0->Set(vtkDataObject::SIL(), this->GetSIL());
  //
  return result;
}
//----------------------------------------------------------------------------
int vtkCircuitReader::OpenReportFile()
{
  std::string reportname = "";
  bbp::Reports_Specification &reports = this->Experiment.reports();
  for (bbp::Reports_Specification::iterator ri=reports.begin(); ri!=reports.end(); ++ri) {
    reportname = (*ri).label();
  }                                          
  if (reports.size()==0) {
    return 0;
  }
  bbp::Reports_Specification::iterator ri=reports.find(reportname);
  this->ReportReader = bbp::CompartmentReportReader::createReader(*ri);
  this->ReportReader->updateMapping(Target);
  //
  this->startTime = (*ri).start_time();
  this->stopTime  = (*ri).end_time();
  this->timestep  = (*ri).timestep();

  // the mapping array(s) provided by the report reader
  this->ReportMapping = this->ReportReader->getMapping();  

  // we need a cell Target object (another neuron list), so get one from the neuron list (waste of memory?)
  bbp::Cell_Target ctarget = this->Target.cell_target();
  if (ctarget.size() != 0) {
    Cell_Index index = 0;
    for (Cell_Target::const_iterator i=ctarget.begin(); i!=ctarget.end(); ++i, ++index) {
      this->OffsetMapping[*i] = index;
    }
  } else {
    bbp::Neurons &neurons = this->Microcircuit->neurons(); 
    Cell_Index j = 0;
    for (bbp::Neurons::const_iterator i = neurons.begin(); i != neurons.end(); ++i,++j) {
      std::cerr << "We must provide some kind of mapping for neuron list " << std::endl;
      //            (*i)->setSimulationBufferIndex(j);
    }
  }

  return 1;
}
//----------------------------------------------------------------------------
// We know that our data contains only triangles, so this is safe
//----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> UnstructuredGridToPolyData(vtkUnstructuredGrid *ug, vtkSmartPointer<vtkPolyData> pd) 
{
  if (!pd) pd = vtkSmartPointer<vtkPolyData>::New();
  pd->SetPoints(ug->GetPoints());
  pd->SetPolys(ug->GetCells());
  pd->GetPointData()->ShallowCopy(ug->GetPointData());
  //
  return pd;
}
//----------------------------------------------------------------------------
int vtkCircuitReader::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  if (!this->FileName) {
    vtkErrorMacro(<< "A BlueConfig FileName must be specified.");
    return 0;
  }

  // get the info object
  vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
//  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

  // get the ouptut
  vtkPointSet *output0 = vtkPointSet::SafeDownCast(outInfo0->Get(vtkDataObject::DATA_OBJECT()));
//  vtkPointSet *output1 = vtkPointSet::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));

  // Which time step has been requested
  double requestedTimeValue = outInfo0->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) 
    ? outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) : 0.0; 
//    : outInfo1->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  //
  this->ActualTimeStep = std::find_if(
    this->TimeStepValues.begin(), this->TimeStepValues.end(),
    std::bind2nd( vtkCircuitReaderToleranceCheck( 
        this->IntegerTimeStepValues ? 0.5 : this->TimeStepTolerance ), requestedTimeValue ))
    - this->TimeStepValues.begin();
  //
  bool NeedToRegernerateTime = false;
  if (requestedTimeValue!=this->CurrentTime) {
    this->CurrentTime = requestedTimeValue;
    NeedToRegernerateTime = true;
  }
  //
  output0->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);
//  output1->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);

  // parallel pieces info
  this->UpdatePiece = outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  vtkSmartPointer<vtkTimerLog> load_timer = vtkSmartPointer<vtkTimerLog>::New();        
  load_timer->StartTimer();
 
  bbp::Neurons &neurons = this->Microcircuit->neurons(); 
  if (neurons.size()==0) {
    // Load neurons for this target so we can partition them
    this->Microcircuit->load(this->Target, bbp::NEURONS);
  }
  int WholeExtent[6] = { 0, neurons.size(), 0, 0, 0, 0 };
  if (this->MaximumNumberOfNeurons>0) {
    WholeExtent[1] = std::min(neurons.size(), (size_t)(this->MaximumNumberOfNeurons));
  }
  vtkSmartPointer<vtkExtentTranslator> extTran = vtkSmartPointer<vtkExtentTranslator>::New();
  extTran->SetSplitModeToBlock();
  extTran->SetNumberOfPieces(this->UpdateNumPieces);
  extTran->SetPiece(this->UpdatePiece);
  extTran->SetWholeExtent(WholeExtent);
  extTran->PieceToExtent();
  extTran->GetExtent(this->PartitionExtents);
  //
  //
  // Create iterators for the begin and end of our partition
  bbp::Neurons::iterator NeuronStart = neurons.begin();
  bbp::Neurons::iterator NeuronEnd   = neurons.begin();
  std::advance(NeuronStart,PartitionExtents[0]);
  std::advance(NeuronEnd,PartitionExtents[1]);
  // create a new target based on our subrange of neurons, clear any contests first.
  this->Partitioned_target = bbp::Target();
  for (bbp::Neurons::iterator neuron=NeuronStart; neuron!=NeuronEnd; ++neuron) {
    this->Partitioned_target.insert(neuron->gid());
  }
  // Load morphology and meshes for this subtarget
#ifdef MANUAL_MESH_LOAD
  this->Microcircuit->load(this->Partitioned_target, bbp::NEURONS | bbp::MORPHOLOGIES); // | bbp::MESHES);
#else
  this->Microcircuit->load(this->Partitioned_target, bbp::NEURONS | bbp::MORPHOLOGIES | bbp::MESHES);
#endif
  //
  bool NeedToRegenerateMesh = (FileModifiedTime>MeshGeneratedTime) || (MeshParamsModifiedTime>MeshGeneratedTime);
  if (NeedToRegenerateMesh) {
    if (this->ExportNeuronMesh) {
      this->GenerateNeuronMesh(request, inputVector, outputVector);
    }
    else {
      this->CachedNeuronMesh->Initialize();
    }
    if (this->ExportMorphologySkeleton) {
      this->GenerateMorphologySkeleton(request, inputVector, outputVector);
    }
    else {
      this->CachedMorphologySkeleton->Initialize();
    }
  
    if (this->UpdateNumPieces>1 && this->ParallelRedistribution) {
      vtkSmartPointer<vtkTimerLog> redist_timer = vtkSmartPointer<vtkTimerLog>::New();        
      redist_timer->StartTimer();
      //
      this->DistributedDataFilter = vtkSmartPointer<vtkDistributedDataFilter>::New();
      this->DistributedDataFilter->SetInputData(this->CachedNeuronMesh);
      this->DistributedDataFilter->SetBoundaryModeToAssignToOneRegion();
      this->DistributedDataFilter->SetClipCells(0);
      this->DistributedDataFilter->SetUseMinimalMemory(1);
      this->DistributedDataFilter->Update();

      // setup bounds for future use
      this->BoundsTranslator->SetKdTree(this->DistributedDataFilter->GetKdtree());

      vtkDataSet *dataset = vtkDataSet::SafeDownCast(this->DistributedDataFilter->GetOutput());
      double bounds[6];
      dataset->GetBounds(bounds);
      this->BoundsTranslator->ExchangeBoundsForAllProcesses(this->Controller, bounds);
      this->BoundsTranslator->InitWholeBounds();
      int whole_extent[6] = {0, 8191, 0, 8191, 0, 8191};
      this->BoundsTranslator->SetWholeExtent(whole_extent);

      this->CachedNeuronMesh = UnstructuredGridToPolyData(vtkUnstructuredGrid::SafeDownCast(this->DistributedDataFilter->GetOutput()), this->CachedNeuronMesh);
      this->DistributedDataFilter->SetInputData(NULL);
      this->DistributedDataFilter = NULL;
      //
      redist_timer->StopTimer();
      if (this->UpdatePiece==0) {
        std::cout << "ParallelRedistribution : " << redist_timer->GetElapsedTime() << " seconds\n";
      }
    }
    this->MeshGeneratedTime.Modified();
  }
  if (NeedToRegenerateMesh || NeedToRegernerateTime) {
    bool do_rep = (this->NumberOfTimeSteps>0) && this->GetPointArrayStatus(BBP_ARRAY_NAME_VOLTAGE);
    if (do_rep && this->UpdateNumPieces==1) {
     this->CreateReportScalars(request, inputVector, outputVector);
    }
    else {
      this->CachedNeuronMesh->GetPointData()->RemoveArray(BBP_ARRAY_NAME_VOLTAGE);
    }
    this->TimeModifiedTime.Modified();
  }
  output0->ShallowCopy(this->CachedNeuronMesh);
//  output1->ShallowCopy(this->CachedMorphologySkeleton);
  //  
  load_timer->StopTimer();
  if (this->UpdatePiece==0) {
    std::cout << "Mesh Load and Redistribution : " << load_timer->GetElapsedTime() << " seconds\n";
  }
  return 1;
}
//-----------------------------------------------------------------------------
void vtkCircuitReader::AddOneNeuronToMesh(bbp::Neuron *neuron, const bbp::Mesh *mesh, vtkIdType Ncount, vtkPoints *points, vtkIdType *cells, vtkFieldData *field, vtkIdType &offsetN, vtkIdType &offsetC)
{  
  bool do_nrm = this->GetPointArrayStatus(BBP_ARRAY_NAME_NORMAL);
  bool do_sid = this->GetPointArrayStatus(BBP_ARRAY_NAME_SECTION_ID);
  bool do_nid = this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONGID);
  bool do_nix = this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONINDEX);
  bool do_rtn = this->GetPointArrayStatus(BBP_ARRAY_NAME_RTNEURON_OPACITY);
  
  //
  vtkSmartPointer<vtkIntArray>   neuronId = do_nid ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_NEURONGID)) : NULL;
  vtkSmartPointer<vtkIntArray>   neuronIx = do_nix ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_NEURONINDEX)) : NULL;
  vtkSmartPointer<vtkIntArray>  sectionId = do_sid ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_SECTION_ID)) : NULL;
  vtkSmartPointer<vtkFloatArray> nvectors = do_nrm ? vtkFloatArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_NORMAL)) : NULL;
  vtkSmartPointer<vtkFloatArray> rtneuron = do_rtn ? vtkFloatArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_RTNEURON_OPACITY)) : NULL;
  //
  bbp::Vertex_Index         vertexCount = mesh->vertex_count();
  bbp::Triangle_Index         faceCount = mesh->triangle_count();
  const Array<Vector_3D<bbp::Micron> >       &vertices = mesh->vertices();
  const Section_ID                        *section_ids = mesh->vertex_sections().pointer();
  const float                               *positions = mesh->vertex_relative_distances().pointer();
  const Array<Vector_3D<bbp::Micron> > &vertex_normals = mesh->normals();
  const bbp::Morphology &morph = neuron->morphology();
    //
  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->PostMultiply();
  transform->RotateY(neuron->orientation().rotation);
  transform->Translate(neuron->position().vector());
  transform->Update();
  vtkIdType insertN = offsetN;
  for (bbp::Vertex_Index v=0 ; v<vertexCount; ++v) {
    float newPoint[3];
    transform->TransformPoint(vertices[v].vector(),newPoint); 
    points->SetPoint(offsetN, newPoint);
    Section_ID sectionID = section_ids[v];

    if (do_nrm) {
      transform->TransformNormal(vertex_normals[v].vector(),newPoint); 
      nvectors->SetTuple(offsetN, newPoint); 
    }
    if (do_nid) {
      neuronId->SetValue(offsetN,neuron->gid());
    }
    if (do_nix) {
      neuronIx->SetValue(offsetN,Ncount);
    }
    if (do_sid) {
      sectionId->SetValue(offsetN, sectionID);
    }
    //
    // Following code taken from RTNeuron
    // C:\Code\Buildyard\src\RTNeuron\src\render\rawModels.cpp
    //
    const Soma_Surface_Points &soma = morph.soma().surface_points();
    if (do_rtn) {
      float position = positions[v];
      const Section &section = morph.section(sectionID);
      if (position < 0) {
              position = 0;
      } else if (position > 1) {
              position = 1;
      }
      float width;
      const float trunkWidth = 2;
      if (section.type() == SOMA) {
          float radius = soma.max_radius() * 0.9;
          float distance = (vertices[v] - soma.center()).length() - radius;
          if (distance < 0) {
              width = radius;
          } else if (distance > 2) {
              width = trunkWidth;
          } else {
              width = ((1 - distance / 2) * radius + 
                        distance / 2 * trunkWidth);
          }
      } else if (section.parent().type() == SOMA) {
          float diameter = section.cross_section(position).diameter();
          float distance = section.length() * position;
          if (distance > 2)
              width = diameter;
          else {
              if (diameter > trunkWidth) {
                  width = diameter / 2;
              } else {
                  width = trunkWidth * (1 - distance / 2) + 
                      diameter * distance / 2;
              }
          }
      } else {
          width = section.cross_section(position).diameter();
      }
      float alpha = 0.35 * (1 - exp(-width * 0.3));
      rtneuron->SetValue(offsetN, alpha);
    }
    offsetN++;
  }

  //
  // Generate triangles in cell array for each face
  // 
  const bbp::Vertex_Index* faces2 = mesh->triangles().pointer();
  for (bbp::Triangle_Index t=0; t<faceCount; ++t)
  {
    cells[offsetC*4 + 0] = 3;
    cells[offsetC*4 + 1] = insertN + faces2[t*3 + 0];
    cells[offsetC*4 + 2] = insertN + faces2[t*3 + 1];
    cells[offsetC*4 + 3] = insertN + faces2[t*3 + 2];
    offsetC++;
  }
}
//-----------------------------------------------------------------------------
void vtkCircuitReader::AddOneNeuronToMorphologySkeleton(bbp::Neuron *neuron, vtkIdType Ncount, vtkPoints *points, vtkIdType *cells, vtkFieldData *field, vtkIdType &offsetN, vtkIdType &offsetC)
{
  bool do_nid = this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONGID); 
  bool do_nix = this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONINDEX);
  bool do_sid = this->GetPointArrayStatus(BBP_ARRAY_NAME_SECTION_ID);
  bool do_sty = this->GetPointArrayStatus(BBP_ARRAY_NAME_SECTION_TYPE);
  bool do_ddr = this->GetPointArrayStatus(BBP_ARRAY_NAME_DENDRITE_RADIUS);
  //
  vtkSmartPointer<vtkIntArray>                neuronId = do_nid ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_NEURONGID)) : NULL;
  vtkSmartPointer<vtkIntArray>                neuronIx = do_nix ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_NEURONINDEX)) : NULL;
  vtkSmartPointer<vtkIntArray>               sectionId = do_sid ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_SECTION_ID)) : NULL;
  vtkSmartPointer<vtkUnsignedCharArray>    sectionType = do_sty ? vtkUnsignedCharArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_SECTION_TYPE)) : NULL;
  vtkSmartPointer<vtkFloatArray>        dendriteRadius = do_ddr ? vtkFloatArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_DENDRITE_RADIUS)) : NULL;
  //
  const bbp::Morphology               *morph = &neuron->morphology();
  Morphology_Dataset                 dataset = /*std::move(*/morph->operator Morphology_Dataset();//);
  const Section_Type          *section_types = dataset.section_types();
  const Transform_3D<bbp::Micron> &transform = neuron->global_transform();

  vtkIdType insertN = offsetN;
  vtkIdType i = 0;
  //
  const bbp::Sections &sections = morph->neurites();
  Segments::const_iterator current_segment;
  for (Sections::const_iterator section=sections.begin(); section!=sections.end(); ++section) {
    Segments segments = section->segments();
    // allowed types are SOMA/AXON/DENDRITE/APICAL_DENDRITE
    bbp::Section_Type section_type = section_types[section->id()];

    // if the segment has more than zero pieces, add the start point
    if (segments.begin()!=segments.end()) {
      bbp::Vector_3D<bbp::Micron> newPoint = transform*segments.begin()->begin().center();
      double radius = segments.begin()->begin().diameter()/2.0; 
      //
      points->SetPoint(offsetN + i, newPoint.vector());
      if (do_sid) sectionId->SetValue(offsetN + i, section->id());
      if (do_sty) sectionType->SetValue(offsetN + i, section_type);
      if (do_ddr) dendriteRadius->SetValue(offsetN + i, radius);
      if (do_nid) neuronId->SetValue(offsetN + i, neuron->gid());
      if (do_nix) neuronIx->SetValue(offsetN,Ncount);
      i++;
    }
    // add one point for each segment piece
    for (Segments::const_iterator segment=segments.begin(); segment!=segments.end(); ++segment) {
      bbp::Vector_3D<bbp::Micron> newPoint = transform*segment->center();
      double radius = segment->diameter()/2.0; 
      //
      points->SetPoint(offsetN + i, newPoint.vector());
      if (do_sid) sectionId->SetValue(offsetN + i, section->id());
      if (do_sty) sectionType->SetValue(offsetN + i, section_type);
      if (do_ddr) dendriteRadius->SetValue(offsetN + i, radius);
      if (do_nid) neuronId->SetValue(offsetN + i, neuron->gid());
      if (do_nix) neuronIx->SetValue(offsetN,Ncount);
      //
      if (i>0) {
        cells[offsetC*3 + 0] = 2;
        cells[offsetC*3 + 1] = i - 1 + insertN;
        cells[offsetC*3 + 2] = i     + insertN;
        offsetC++;
      }
      i++;
    }
  }
  offsetN += i;
}
//-----------------------------------------------------------------------------
void vtkCircuitReader::GenerateNeuronMesh(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
//  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

  // VTK arrays 
  vtkSmartPointer<vtkPoints>           points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray>     triangles = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPointData>     pointdata = vtkSmartPointer<vtkPointData>::New();
  vtkSmartPointer<vtkIntArray>       neuronId;
  vtkSmartPointer<vtkIntArray>       neuronIx;
  vtkSmartPointer<vtkIntArray>      sectionId;
  vtkSmartPointer<vtkFloatArray>     nvectors;
  vtkSmartPointer<vtkFloatArray>     rtneuron;

  // counters
  vtkIdType maxPoints = 0;
  vtkIdType maxCells  = 0;

  //
  // loop over neurons : count up the total vertices and cells so we can allocate memory all in one go
  //
  bbp::Neurons &neurons = this->Microcircuit->neurons(); 

#ifdef MANUAL_MESH_LOAD
  //
  bbp::URI location_m = this->Microcircuit->reader().data_source("mesh");
  bbp::URI location_c = this->Microcircuit->reader().data_source("circuit");
  Mesh_Reader_Ptr mesh_reader = Mesh_Reader::create_reader(location_m);
  Microcircuit_Composition_Reader_Ptr circuit_reader = Microcircuit_Composition_Reader::create_reader(location_c);
  typedef boost::shared_ptr<bbp::Meshes> Meshes_Ptr;

  bbp::Meshes meshes;
  // Loading meshes
  mesh_reader->load(
      meshes, this->Partitioned_target,
      circuit_reader->source(),
      true, // vertices
      true, // triangles
      true, // mapping
      false); // strips
  
  // free unwanted memory
  mesh_reader.reset();
  circuit_reader.reset();
#endif

  for (bbp::Target::cell_iterator gid=this->Partitioned_target.cell_begin(); gid!=this->Partitioned_target.cell_end(); ++gid) {
    Neurons::iterator neuron = neurons.find( *gid );
#ifdef MANUAL_MESH_LOAD
    const bbp::Morphology *morphology = &(*neuron).morphology();
    Meshes::iterator mesh = meshes.find(morphology->label());
    if (mesh != meshes.end()) {
#else
    const bbp::Mesh *mesh = &neuron->morphology().mesh(); 
    {
#endif
      maxPoints += mesh->vertex_count();
      maxCells += mesh->triangle_count();
    }
  }
  //
  // reserve space for coordinates
  //
  points->GetData()->Resize(maxPoints);
  points->SetNumberOfPoints(maxPoints);
  //
  // reserve space for each scalar/vector field
  //
  bool do_nid = this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONGID);
  bool do_nix = this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONINDEX);
  bool do_sid = this->GetPointArrayStatus(BBP_ARRAY_NAME_SECTION_ID);
  bool do_nrm = this->GetPointArrayStatus(BBP_ARRAY_NAME_NORMAL);
  bool do_rtn = this->GetPointArrayStatus(BBP_ARRAY_NAME_RTNEURON_OPACITY);
  //
  if (do_nid) {
    neuronId = vtkSmartPointer<vtkIntArray>::New();
    neuronId->SetName(BBP_ARRAY_NAME_NEURONGID);
    neuronId->Resize(maxPoints);
    neuronId->SetNumberOfTuples(maxPoints);
    pointdata->AddArray(neuronId);
  }
  //
  if (do_nix) {
    neuronIx = vtkSmartPointer<vtkIntArray>::New();
    neuronIx->SetName(BBP_ARRAY_NAME_NEURONINDEX);
    neuronIx->Resize(maxPoints);
    neuronIx->SetNumberOfTuples(maxPoints);
    pointdata->AddArray(neuronIx);
  }
  //    
  if (do_sid) {
    sectionId = vtkSmartPointer<vtkIntArray>::New();
    sectionId->SetName(BBP_ARRAY_NAME_SECTION_ID);
    sectionId->Resize(maxPoints);
    sectionId->SetNumberOfTuples(maxPoints);
    pointdata->AddArray(sectionId);
  }
  //
  if (do_nrm) {
    nvectors = vtkSmartPointer<vtkFloatArray>::New();
    nvectors->SetNumberOfComponents(3);
    nvectors->SetName(BBP_ARRAY_NAME_NORMAL);
    nvectors->Resize(maxPoints);
    nvectors->SetNumberOfTuples(maxPoints);
    pointdata->SetNormals(nvectors);
  }
  //
  if (do_rtn) {
    rtneuron = vtkSmartPointer<vtkFloatArray>::New();
    rtneuron->SetName(BBP_ARRAY_NAME_RTNEURON_OPACITY);
    rtneuron->Resize(maxPoints);
    rtneuron->SetNumberOfTuples(maxPoints);
    pointdata->AddArray(rtneuron);
  }
  //
  // reserve space for cells
  //
  vtkIdType *cells = triangles->WritePointer(maxCells, 4*(maxCells));
  //
  // Each neuron counts vertices starting from zero, but we increment one much larger array
  // Track insertion location for Cell vertex index Ids
  vtkIdType offsetN = 0;
  vtkIdType offsetC = 0;
  vtkIdType Ncount = this->PartitionExtents[0];
  for (bbp::Target::cell_iterator gid=this->Partitioned_target.cell_begin(); gid!=this->Partitioned_target.cell_end(); ++gid, ++Ncount) {
    Neurons::iterator neuron = neurons.find( *gid );
    //
#ifdef MANUAL_MESH_LOAD
    const bbp::Morphology *morphology = &(*neuron).morphology();
    Meshes::iterator mesh = meshes.find(morphology->label());
    if (mesh != meshes.end()) {
#else
    const bbp::Mesh *mesh = &neuron->morphology().mesh(); 
    {
#endif
      this->AddOneNeuronToMesh(&*neuron, &*mesh, Ncount, points, cells, pointdata, offsetN, offsetC);
    }
  }

#ifdef MANUAL_MESH_LOAD
  meshes.clear();
#endif
  // get the outputs
  vtkPolyData *output0 = vtkPolyData::SafeDownCast(outInfo0->Get(vtkDataObject::DATA_OBJECT()));
//  vtkPolyData *output1 = vtkPolyData::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));


  //
  this->CachedNeuronMesh->SetPoints(points);
  this->CachedNeuronMesh->SetPolys(triangles);
  this->CachedNeuronMesh->GetPointData()->ShallowCopy(pointdata);
}
//-----------------------------------------------------------------------------
void vtkCircuitReader::GenerateMorphologySkeleton(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
//  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

  // get the outputs
  vtkPolyData *output0 = vtkPolyData::SafeDownCast(outInfo0->Get(vtkDataObject::DATA_OBJECT()));
//  vtkPolyData *output1 = vtkPolyData::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));

  // VTK arrays 
  vtkSmartPointer<vtkPoints>                    points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray>                  lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPointData>              pointdata = vtkSmartPointer<vtkPointData>::New();
  vtkSmartPointer<vtkFloatArray>        dendriteRadius;
  vtkSmartPointer<vtkUnsignedCharArray>    sectionType;
  vtkSmartPointer<vtkIntArray>               sectionId;
  vtkSmartPointer<vtkIntArray>                neuronId;
  vtkSmartPointer<vtkIntArray>                neuronIx;

  // counters
  vtkIdType maxPoints = 0;
  vtkIdType maxCells  = 0;

  bbp::Neurons &neurons = this->Microcircuit->neurons(); 
  //
  // loop over neurons : count up the total vertices and cells so we can allocate memory all in one go
  //
  for (bbp::Target::cell_iterator gid=this->Partitioned_target.cell_begin(); gid!=this->Partitioned_target.cell_end(); ++gid) {
    Neurons::iterator neuron = neurons.find( *gid );
    //
    const bbp::Sections &sections = neuron->morphology().neurites();
    for (Sections::const_iterator section = sections.begin(); section != sections.end(); ++section) {
      Segments segments = section->segments();
      if (segments.begin()!=segments.end()) {
        ++maxPoints;
      }
      for (Segments::const_iterator segment = segments.begin(); segment != segments.end(); ++segment) {
        if (maxPoints>0) {
          maxCells++;
        }
        ++maxPoints;
      }
    }
  }
  //
  // reserve space for coordinates
  //
  points->GetData()->Resize(maxPoints);
  points->SetNumberOfPoints(maxPoints);
  //
  // reserve space for each scalar/vector field
  //
  bool do_ddr = this->GetPointArrayStatus(BBP_ARRAY_NAME_DENDRITE_RADIUS);
  bool do_sty = this->GetPointArrayStatus(BBP_ARRAY_NAME_SECTION_TYPE);
  bool do_sid = this->GetPointArrayStatus(BBP_ARRAY_NAME_SECTION_ID);
  bool do_nid = this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONGID); 
  bool do_nix = this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONINDEX);
  //
  if (do_nid) {
    neuronId = vtkSmartPointer<vtkIntArray>::New();
    neuronId->SetName(BBP_ARRAY_NAME_NEURONGID);
    neuronId->Resize(maxPoints);
    neuronId->SetNumberOfTuples(maxPoints);
    pointdata->AddArray(neuronId);
  }
  //
  if (do_nix) {
    neuronIx = vtkSmartPointer<vtkIntArray>::New();
    neuronIx->SetName(BBP_ARRAY_NAME_NEURONINDEX);
    neuronIx->Resize(maxPoints);
    neuronIx->SetNumberOfTuples(maxPoints);
    pointdata->AddArray(neuronIx);
  }
  //    
  if (do_sid) {
    sectionId = vtkSmartPointer<vtkIntArray>::New();
    sectionId->SetName(BBP_ARRAY_NAME_SECTION_ID);
    sectionId->Resize(maxPoints);
    sectionId->SetNumberOfTuples(maxPoints);
    pointdata->AddArray(sectionId);
  }
  //
  if (do_sty) {
    sectionType = vtkSmartPointer<vtkUnsignedCharArray>::New();
    sectionType->SetName(BBP_ARRAY_NAME_SECTION_TYPE);
    sectionType->Resize(maxPoints);
    sectionType->SetNumberOfTuples(maxPoints);
    pointdata->AddArray(sectionType);
  }
  //
  if (do_ddr) {
    dendriteRadius = vtkSmartPointer<vtkFloatArray>::New();
    dendriteRadius->SetName(BBP_ARRAY_NAME_DENDRITE_RADIUS);
    dendriteRadius->Resize(maxPoints);
    dendriteRadius->SetNumberOfTuples(maxPoints);
    pointdata->AddArray(dendriteRadius);
  }
  //
  // reserve space for cells
  //
  vtkIdType *cells = lines->WritePointer(maxCells, 3*(maxCells));
  //
  // Each neuron counts vertices starting from zero, but we increment one much larger array
  // Track insertion location for Cell vertex index Ids
  vtkIdType offsetN = 0;
  vtkIdType offsetC = 0;

  vtkIdType Ncount = this->PartitionExtents[0];
  for (bbp::Target::cell_iterator gid=this->Partitioned_target.cell_begin(); gid!=this->Partitioned_target.cell_end(); ++gid, ++Ncount) {
    Neurons::iterator neuron = neurons.find( *gid );
    //
    this->AddOneNeuronToMorphologySkeleton(&*neuron, Ncount, points, cells, pointdata, offsetN, offsetC);
  }
  //
  this->CachedMorphologySkeleton->SetPoints(points);
  this->CachedMorphologySkeleton->SetLines(lines);
  this->CachedMorphologySkeleton->GetPointData()->ShallowCopy(pointdata);
}
//-----------------------------------------------------------------------------
void vtkCircuitReader::CreateReportScalars(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  this->ReportReader->loadFrame(this->CurrentTime, this->_currentFrame);
  // bind the current frame to the microcircuit
  this->Microcircuit->update( this->_currentFrame );

  // Voltage level provided by simulation reports (Morphology Skeleton)
  vtkSmartPointer<vtkFloatArray> voltageM;
  if (this->ExportMorphologySkeleton) {
    voltageM = vtkSmartPointer<vtkFloatArray>::New();
    voltageM->SetName(BBP_ARRAY_NAME_VOLTAGE);
    this->CachedMorphologySkeleton->GetPointData()->AddArray(voltageM);
    vtkIdType maxPoints = this->CachedMorphologySkeleton->GetPoints()->GetNumberOfPoints();
    voltageM->Resize(maxPoints);
    voltageM->SetNumberOfTuples(maxPoints);
  }

  // Voltage level provided by simulation reports (Neuron Mesh)
  vtkSmartPointer<vtkFloatArray> voltageN = vtkSmartPointer<vtkFloatArray>::New();
  if (this->ExportNeuronMesh) {
    voltageN = vtkSmartPointer<vtkFloatArray>::New();
    voltageN->SetName(BBP_ARRAY_NAME_VOLTAGE);
    this->CachedNeuronMesh->GetPointData()->AddArray(voltageN);
    vtkIdType maxPointsN = this->CachedNeuronMesh->GetPoints()->GetNumberOfPoints();
    voltageN->Resize(maxPointsN);
    voltageN->SetNumberOfTuples(maxPointsN);
  }

  bbp::Neurons &neurons = this->Microcircuit->neurons(); 
  vtkIdType offsetN = 0, offsetM = 0;
  for (bbp::Target::cell_iterator gid=this->Partitioned_target.cell_begin(); gid!=this->Partitioned_target.cell_end(); ++gid) {
    Neurons::iterator neuron = neurons.find( *gid );
    //
    if (this->ExportMorphologySkeleton) {
      offsetM = this->AddReportScalarsToNeuronMorphology(&*neuron, voltageM, offsetM);  
    }
    if (this->ExportNeuronMesh) {
      offsetN = this->AddReportScalarsToNeuronMesh(&*neuron, voltageN, offsetN);  
    }
  }
}
//-----------------------------------------------------------------------------
vtkIdType vtkCircuitReader::AddReportScalarsToNeuronMorphology(bbp::Neuron *neuron, vtkFloatArray *voltages, vtkIdType offsetN)
{
  const bbp::Morphology* morph                    = &neuron->morphology();
  const bbp::Transform_3D<bbp::Micron> &transform = neuron->global_transform();
  const bbp::Count index = this->OffsetMapping[neuron->gid()];
  std::cout << "Neuron with GID " << neuron->gid() << " has mapping offset " << index << std::endl;

  Morphology_Dataset dataset = /*std::move(*/morph->operator Morphology_Dataset();//);
  const Section_Type                *section_types = dataset.section_types();

  /* These are the buffer offsets for this neuron */
  const std::vector<Report_Frame_Index> &offsets = this->ReportMapping->sections_offsets(index);

  const float *buffer = this->_currentFrame.getData<bbp::Voltages>().get();
  float somaVoltage = buffer[offsets[0]];
  float rvoltage = -70;

  // Finding the voltage value of the last compartment of the last axon section. 
  // This is not completely accurate but it's the best we can do to assign some color 
  // to axon sections with undefined simulation data.
  const Sections &axon = morph->axon();
  const Section_ID lastAxon = axon.size() == 0 || this->ReportMapping->number_of_compartments(index, axon.begin()->id()) == 0 ? axon.begin()->id() : morph->soma().id();
  float undefinedAxonVoltage = buffer[offsets[lastAxon] + this->ReportMapping->number_of_compartments(index, lastAxon)];

  vtkIdType i = 0;
  const bbp::Sections &sections = morph->neurites();
  Segments::const_iterator current_segment;
  for (Sections::const_iterator section=sections.begin(); section!=sections.end(); ++section) {
    Segments segments = section->segments();
    // allowed types are SOMA/AXON/DENDRITE/APICAL_DENDRITE
    bbp::Section_Type section_type = section_types[section->id()];

    // if the segment has more than zero pieces, add the start point
    if (segments.begin()!=segments.end()) {
      bbp::Vector_3D<bbp::Micron> newPoint = transform*segments.begin()->begin().center();
      double radius = segments.begin()->begin().diameter()/2.0; 
      //
      uint16_t compartments = this->ReportMapping->number_of_compartments(index, section->id());

      if (compartments) {
        // Computing the relative length within section of the capsule midpoint
        float position = section->section_distance(*segments.begin());
        uint16_t compartment = std::min(compartments - 1, (int)floor(compartments * position));
        unsigned long offset = offsets[section->id()];
        rvoltage = buffer[offset + compartment];
      } else {
        rvoltage = -100;
/*
        if (section->type() == SOMA) {
          rvoltage = somaVoltage;
        } else if (section->type() == AXON) {
          rvoltage = undefinedAxonVoltage;
        }
*/
      }
      voltages->SetValue(offsetN + i, rvoltage);
      i++;
    }
    // add one point for each segment piece
    for (Segments::const_iterator segment=segments.begin(); segment!=segments.end(); ++segment) {
      uint16_t compartments = this->ReportMapping->number_of_compartments(index, section->id());
      bbp::Vector_3D<bbp::Micron> newPoint = transform*segment->center();
      double radius = segment->diameter()/2.0; 

      if (compartments) {
        // Computing the relative length within section of the capsule midpoint
        float position = section->section_distance(*segment);
        uint16_t compartment = std::min(compartments - 1, (int)floor(compartments * position));
        unsigned long offset = offsets[section->id()];
        rvoltage = buffer[offset + compartment];
      } else {
        rvoltage = -100;
/*
        if (section->type() == SOMA) {
          rvoltage = somaVoltage;
        } else if (section->type() == AXON) {
          rvoltage = undefinedAxonVoltage;
        }
*/
      }

      voltages->SetValue(offsetN + i, rvoltage);
      i++;
    }
  }
  return offsetN + i;
}
//-----------------------------------------------------------------------------
vtkIdType vtkCircuitReader::AddReportScalarsToNeuronMesh(bbp::Neuron *neuron, vtkFloatArray *voltages, vtkIdType offsetN)
{
  const bbp::Morphology *sdk_morph    = &neuron->morphology();
  const bbp::Mesh       *sdk_mesh     = &sdk_morph->mesh();
  //
  const Array<Vector_3D<bbp::Micron> >      &vertices2 = sdk_mesh->vertices();
  const Array<Section_ID>                 &section_ids = sdk_mesh->vertex_sections();
  const Array<float>                &section_distances = sdk_mesh->vertex_relative_distances();
  const Array<Vector_3D<bbp::Micron> > &vertex_normals = sdk_mesh->normals();
  //

  vtkIdType vertexCount = sdk_mesh->vertex_count();

  // get the index 
  Cell_Index cell_index = neuron->index();
  if (cell_index==UNDEFINED_CELL_INDEX) {
    for (bbp::Vertex_Index i=0; i<vertexCount; ++i) {
      voltages->SetValue(offsetN + i, -100);
    }
    return offsetN + vertexCount;
  }

  const bbp::Count index = this->OffsetMapping[neuron->gid()];

  // These are the buffer offsets for this neuron 
  const std::vector<Report_Frame_Index> &offsets = this->ReportMapping->sections_offsets(index);

  const float *buffer = this->_currentFrame.getData<bbp::Voltages>().get();
  float somaVoltage = buffer[offsets[0]];
  float rvoltage = -70;

  // Finding the voltage value of the last compartment of the last axon section. 
  // This is not completely accurate but it's the best we can do to assign some color 
  // to axon sections with undefined simulation data.
  const Sections &axon = sdk_morph->axon();
  const Section_ID lastAxon = axon.size() == 0 || this->ReportMapping->number_of_compartments(cell_index, axon.begin()->id()) == 0 ? axon.begin()->id() : sdk_morph->soma().id();
  float undefinedAxonVoltage = buffer[offsets[lastAxon] + this->ReportMapping->number_of_compartments(cell_index, lastAxon)];

  for (bbp::Vertex_Index i=0; i<vertexCount; ++i) {
    Section_ID     sectionId = section_ids[i];
    size_t num_sections = this->ReportMapping->number_of_sections(index);
    uint16_t compartments = 0;
    if (sectionId < num_sections) {
      compartments = this->ReportMapping->number_of_compartments(index, sectionId);
    }
    float position = section_distances[i];
    //
    if (compartments) {
      // Computing the relative length within section of the capsule midpoint
      uint16_t compartment = std::min(compartments - 1, (int)floor(compartments * position));
      unsigned long offset = offsets[sectionId];
      rvoltage = buffer[offset + compartment];
    } else {
        rvoltage = -100;
/*
        if (section->type() == SOMA) {
          rvoltage = somaVoltage;
        } else if (section->type() == AXON) {
          rvoltage = undefinedAxonVoltage;
        }
*/
    }
    voltages->SetValue(offsetN + i, rvoltage);
  }
  return offsetN + vertexCount;
}
//-----------------------------------------------------------------------------
void vtkCircuitReader::BuildSIL()
{
  // Initialize the SIL, dump all previous information.
  this->SIL->Initialize();

  vtkSmartPointer<vtkVariantArray> childEdge = vtkSmartPointer<vtkVariantArray>::New();
  childEdge->InsertNextValue(0);

  vtkSmartPointer<vtkVariantArray> crossEdge = vtkSmartPointer<vtkVariantArray>::New();
  crossEdge->InsertNextValue(0);

  // CrossEdge is an edge linking hierarchies.
  vtkSmartPointer<vtkUnsignedCharArray> crossEdgesArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  crossEdgesArray->SetName("CrossEdges");
  this->SIL->GetEdgeData()->AddArray(crossEdgesArray);

  std::deque<std::string> names;
  std::deque<vtkIdType> crossEdgeIds; // ids for edges between trees.
  int cc;

  // Now build the hierarchy.
  vtkIdType rootId = this->SIL->AddVertex();
  names.push_back("SIL");

  // Add the Heirarchy subtree.
  vtkIdType targetsRoot = this->SIL->AddChild(rootId, childEdge);
  names.push_back("Targets");

  // Get default targets for the microcircuit.
  if (this->Microcircuit) {
    const bbp::Targets &default_targets = this->Experiment.targets();
    for (bbp::Targets::const_iterator ti=default_targets.begin(); ti!=default_targets.end(); ++ti) {
      std::string name = (*ti).name();
      vtkIdType childBlock = this->SIL->AddChild(targetsRoot, childEdge);
      names.push_back(name.c_str());
      this->TargetsSelection->AddArray(name.c_str());
    }

    // Get user targets for the microcircuit.
    const bbp::Targets &user_targets = this->Experiment.user_targets();
    for (bbp::Targets::const_iterator ti=user_targets.begin(); ti!=user_targets.end(); ++ti) {
      std::string name = (*ti).name();
      vtkIdType childBlock = this->SIL->AddChild(targetsRoot, childEdge);
      names.push_back(name.c_str());
      this->TargetsSelection->AddArray(name.c_str());
    }
  }

  // This array is used to assign names to nodes.
  vtkSmartPointer<vtkStringArray> namesArray = vtkSmartPointer<vtkStringArray>::New();
  namesArray->SetName("Names");
  namesArray->SetNumberOfTuples(this->SIL->GetNumberOfVertices());
  this->SIL->GetVertexData()->AddArray(namesArray);

  std::deque<std::string>::iterator iter;
  for (cc=0, iter = names.begin(); iter != names.end(); ++iter, ++cc)
  {
    const char *name = (*iter).c_str();
    namesArray->SetValue(cc, name);
  }
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
//-----------------------------------------------------------------------------
vtkSmartPointer<vtkMutableDirectedGraph> vtkCircuitReader::GetSIL()
{
  return this->SIL;
}
//----------------------------------------------------------------------------
int vtkCircuitReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
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
    this->MeshParamsModifiedTime.Modified();
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
int vtkCircuitReader::GetNumberOfTargets()
{
  return this->TargetsSelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
const char* vtkCircuitReader::GetTargetsName(int index)
{
  return this->TargetsSelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkCircuitReader::GetTargetsStatus(const char* name)
{
  return this->TargetsSelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkCircuitReader::SetTargetsStatus(const char* name, int status)
{
  if(status)
  {
    this->TargetsSelection->EnableArray(name);
  }
  else
  {
    this->TargetsSelection->DisableArray(name);
  }
}
//----------------------------------------------------------------------------
void vtkCircuitReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
    << (this->FileName ? this->FileName : "(none)") << "\n";  
}
//----------------------------------------------------------------------------
