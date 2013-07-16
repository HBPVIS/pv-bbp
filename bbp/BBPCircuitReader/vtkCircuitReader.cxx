// Visual studio debug settings for paths specific to this module
//
// PATH=D:\build\paraview-3.98\bin\Debug;C:\Program Files\hdf5-vfd-1.8.9\bin;%PATH%
// PV_PLUGIN_PATH=D:\build\buildyard\ParaBBP\bin\Debug
// _NO_DEBUG_HEAP=1
// working directory : D:\build\paraview-3.98\bin\Debug
//
//
// For PARAVIEW_USE_MPI
#include "vtkPVConfig.h"
#ifdef PARAVIEW_USE_MPI
#include "vtkMPI.h"
#include "vtkMPIController.h"
#include "vtkMPICommunicator.h"
#endif
#include "vtkDummyController.h"
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
//
#include "vtkPKdTree.h"
#ifdef PV_BBP_USE_ZOLTAN
 #include "vtkBoundsExtentTranslator.h"
 #include "vtkMeshPartitionFilter.h"
#endif
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
//#include "BBP/IO/File/Parsers/Morphology_HDF5_File_Parser.h"
//#include "BBP/Microcircuit/Readers/Microcircuit_Reader.h"
//#include "BBP/Microcircuit/Readers/Mesh_Reader.h"
//#include "BBP/Microcircuit/Readers/compartmentReportReader.h"
//#include "BBP/Microcircuit/Mappings/Compartment_Report_Mapping.h"
//#include "BBP/Microcircuit/Soma.h"
//#include "BBP/Microcircuit/Mesh.h"
//#include "BBP/Microcircuit/Datasets/Morphology_Dataset.h"

// Header of this Reader
#include "vtkCircuitReader.h"
//#include "SpikeData.h"

// BBP-SDK
// Voxelization
//#include "BBP/Voxelization/voxelization.h"
//#include "BBP/VtkDebugging/visualization.h"
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#define BOOL(x) (x!=0)
//----------------------------------------------------------------------------
#define JB_DEBUG__
#ifdef JB_DEBUG__
#undef  OUTPUTTEXT
#define OUTPUTTEXT(a) std::cout << (a) << std::endl; std::cout.flush();

#undef vtkDebugMacro
#define vtkDebugMacro(a)  \
    { \
    if (this->UpdatePiece>=0) { \
    vtkOStreamWrapper::EndlType endl; \
    vtkOStreamWrapper::UseEndl(endl); \
    vtkOStrStreamWrapper vtkmsg; \
    vtkmsg << "P(" << this->UpdatePiece << "): " a << "\n"; \
    OUTPUTTEXT(vtkmsg.str()); \
    vtkmsg.rdbuf()->freeze(0); \
    } \
    }

#undef vtkErrorMacro
#define vtkErrorMacro(a) vtkDebugMacro(a)  
#endif
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//#define MANUAL_MESH_LOAD
//----------------------------------------------------------------------------
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
vtkCircuitReader::vtkCircuitReader() :
  PrimaryTarget("dummy", bbp::TARGET_CELL), 
  Partitioned_target("dummy", bbp::TARGET_CELL)
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
  this->DeleteExperiment                = 1;
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
#ifdef PV_BBP_USE_ZOLTAN
  this->MeshPartitionFilter           = NULL;
  this->BoundsTranslator              = vtkSmartPointer<vtkBoundsExtentTranslator>::New();
#endif
  //
  this->HyperPolarizedVoltage    = -85.0;
  this->DePolarizedVoltage       = -50.0;
  this->RestingPotentialVoltage  = -65.0;

  this->NumberOfPointsBeforePartitioning = 0;
}
//----------------------------------------------------------------------------
vtkCircuitReader::~vtkCircuitReader()
{
  this->SIL                           = NULL;
  this->CachedNeuronMesh              = NULL;
  this->CachedMorphologySkeleton      = NULL;
#ifdef PV_BBP_USE_ZOLTAN
  this->MeshPartitionFilter           = NULL;
  this->BoundsTranslator              = NULL;
#endif
  //
  this->PointDataArraySelection->FastDelete();
  this->TargetsSelection->FastDelete();
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
int vtkCircuitReader::RequestTimeInformation(
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
}
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
  bool NeedToReloadFile = (FileModifiedTime>FileOpenedTime);

  if (!vtksys::SystemTools::FileExists(this->FileName)) {
    vtkWarningMacro("File not found " << this->FileName);
    NeedToReloadFile = 0;
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

  bool ok = true;
  if (NeedToReloadFile) {
    std::string blueconfig = this->FileName;
    // -------------------------------------------------------------------   
    // Create BBP-SDK Experiment and Microcircuit to access to the neurons.
    // -------------------------------------------------------------------   
    try {
      this->Experiment.close();
      this->Experiment.open(blueconfig);
      this->Microcircuit = this->Experiment.microcircuit_ptr();
    }
    catch (...) {
      this->Experiment.clear();
      std::cout <<"An exception occurred opening the bbp::Experiment " << std::endl;
    }
    //
    this->BuildSIL();
    outInfo0->Set(vtkDataObject::SIL(), this->GetSIL());
    this->FileOpenedTime.Modified();
  }


  ok = true;
  bool NeedToRegenerateInfo = /*NeedToReloadFile || */(TargetsModifiedTime>InfoGeneratedTime);
  if (NeedToRegenerateInfo) {
    // default Target?
    this->TargetName = this->DefaultTarget ? this->DefaultTarget : this->TargetsSelection->GetArrayName(0);
    //
    this->PrimaryTarget = bbp::Target("empty",bbp::TARGET_CELL);
    //
    try {
      int N = this->TargetsSelection->GetNumberOfArrays();
      for (int i=0; i<N; i++) {
        const char *name = this->TargetsSelection->GetArrayName(i);
        if (this->TargetsSelection->ArrayIsEnabled(name)) {
          //          this->TargetName = name;
          try {
            try {
              bbp::Target temp = this->Experiment.user_targets().get_target(name);
              if (temp.size()>0) {
                std::cout << "Adding target " << name << " to list" << std::endl;
                this->PrimaryTarget.insert(temp);
              }
            }
            catch (...) {
              bbp::Target temp = this->Experiment.targets().get_target(name);
              if (temp.size()>0) {
                std::cout << "Adding target " << name << " to list" << std::endl;
                this->PrimaryTarget.insert(temp);
              }
            }
          }
          catch (...) {
            std::cout << "Could not add target " << name << " to list" << std::endl;
          }
        }
      }

      // Don't load meshes yet, we'll do that once we've decided which neurons this node will generate
//      this->Microcircuit->load(this->PrimaryTarget, 0); // bbp::NEURONS);
//      bbp::Cell_Target cellTarget = this->PrimaryTarget.cell_target();
    }
    catch (...) {
      vtkErrorMacro(<<"Could not open the circuit");
      this->PrimaryTarget = bbp::Target("exception",bbp::TARGET_CELL);; 
      ok = false;
    }
    this->InfoGeneratedTime.Modified();
  }

  //    std::cout << cellTarget << std::endl;

  // std::cout <<"Made it past load " << std::endl;
  bool needToRegenerateTimeInfo = true;
  if (needToRegenerateTimeInfo) {
    // time steps of reports are in the report file
    if (ok && this->UpdateNumPieces==1) {
      if (this->OpenReportFile()) {
        this->NumberOfTimeSteps = (this->stopTime-this->startTime)/this->timestep;
      }
      else {
        this->NumberOfTimeSteps = 0;
      }
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

#ifdef PV_BBP_USE_ZOLTAN
  if (this->Controller->GetNumberOfProcesses()>1 && this->ParallelRedistribution) {
    outInfo0->Set(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR(), this->BoundsTranslator);
  }
#endif

  //
  return result;
}
//----------------------------------------------------------------------------
int vtkCircuitReader::OpenReportFile()
{
  std::string reportname = "";
  std::string ideal = "voltage1K";
  bbp::Reports_Specification &reports = this->Experiment.reports();
  for (bbp::Reports_Specification::iterator ri=reports.begin(); ri!=reports.end(); ++ri) {
    reportname = (*ri).label();
    if (ideal==reportname) {
      break;
    }
  }                                          
  if (reports.size()==0) {
    return 0;
  }
  bbp::Reports_Specification::iterator ri=reports.find(reportname);
  this->ReportReader.reset(new bbp::CompartmentReportReader(*ri, this->Partitioned_target));


//  this->ReportReader->getCellTarget();
  //
  this->startTime = (*ri).start_time();
  this->stopTime  = (*ri).end_time();
  this->timestep  = (*ri).timestep();

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

  bool NeedToRegenerateMesh = (FileModifiedTime>MeshGeneratedTime) || (MeshParamsModifiedTime>MeshGeneratedTime)
    || (TargetsModifiedTime>MeshGeneratedTime);
  //
  if (NeedToRegenerateMesh) {
    //
//    this->Microcircuit->load(this->PrimaryTarget, 0); // bbp::NEURONS);
//    bbp::Cell_Target cellTarget = this->PrimaryTarget.cell_target();
//    if (neurons.size()==0) {
      // Load neurons for this target so we can partition them
    this->Microcircuit->load(this->PrimaryTarget, bbp::NEURONS);
    bbp::Neurons &neurons = this->Microcircuit->neurons(); 
//    }
    int WholeExtent[6] = { 0, (int)neurons.size(), 0, 0, 0, 0 };
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

    // neurons are ordered in layers and higher layers have bigger cell counts
    // so read them using a random shuffle to avoid one process getting all the 
    // small ones and another the big ones. We can use an operator[]
    // to get neurons from the container because it uses a map, with the GID as key
    // so build a list of all keys and then shuffle that and let each process use
    // a chunk of the list
    // NB. we want all processes to have the same sequence, seed fixed num
    std::srand(12345); 
    std::vector<uint32_t> shufflevector;
    shufflevector.reserve(neurons.size());
    //
    for (bbp::Neurons::iterator neuron=neurons.begin(); neuron!=neurons.end(); ++neuron) {
      shufflevector.push_back(neuron->gid());
    }
    std::random_shuffle ( shufflevector.begin(), shufflevector.end() );
    //  std::ostream_iterator<vtkIdType> out_it(cout,", ");
    //  std::copy(shufflevector.begin(), shufflevector.end(), out_it );

    // create a new target based on our subrange of neurons, clear any contests first.
    this->Partitioned_target = bbp::Target("ParaViewCells", bbp::TARGET_CELL);

    for (int i=PartitionExtents[0]; i<PartitionExtents[1]; i++) {
      uint32_t gid = shufflevector[i];
      Neurons::iterator ni = neurons.find( gid );
      this->Partitioned_target.insert(gid);
      //    vtkDebugMacro(<< "Adding neuron with GID " << gid);
      //    Neurons::iterator ni = neurons.find( gid );
      //    size_t cell_index = ni->index();
      //    if (cell_index==UNDEFINED_CELL_INDEX) {
      //    }
    }

    // Load morphology and meshes for this subtarget
#ifdef MANUAL_MESH_LOAD
    this->Microcircuit->load(this->Partitioned_target, bbp::NEURONS | bbp::MORPHOLOGIES); // | bbp::MESHES);
#else
    this->Microcircuit->load(this->Partitioned_target, bbp::NEURONS | bbp::MORPHOLOGIES | bbp::MESHES);
#endif
    //
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

    this->NumberOfPointsBeforePartitioning = this->CachedNeuronMesh->GetNumberOfPoints();
#ifdef PV_BBP_USE_ZOLTAN
    if (this->UpdateNumPieces>1 && this->ParallelRedistribution) {
      vtkSmartPointer<vtkTimerLog> redist_timer = vtkSmartPointer<vtkTimerLog>::New();        
      redist_timer->StartTimer();
      //
      this->MeshPartitionFilter = vtkSmartPointer<vtkMeshPartitionFilter>::New();
      this->MeshPartitionFilter->SetInputData(this->CachedNeuronMesh);
      // thell the partition filter we can dump the input memory when needed
      this->MeshPartitionFilter->SetInputDisposable(1);
      // for animation over time, keep the map of point send/receive
      // @TODO : Only save this if voltage reports are being loaded
      this->MeshPartitionFilter->SetKeepInversePointLists(1);

      // release our reference count for now
      this->CachedNeuronMesh = NULL;
      // Update the partition filter with parallel information
      vtkStreamingDemandDrivenPipeline::SafeDownCast(this->MeshPartitionFilter->GetExecutive())
        ->SetUpdateExtent(0, this->UpdatePiece, this->UpdateNumPieces, 0);
      this->MeshPartitionFilter->Update();

      // setup bounds for future use
      this->BoundsTranslator->SetKdTree(this->MeshPartitionFilter->GetKdtree());

      vtkDataSet *dataset = vtkDataSet::SafeDownCast(this->MeshPartitionFilter->GetOutput());
      double bounds[6];
      dataset->GetBounds(bounds);
      this->BoundsTranslator->ExchangeBoundsForAllProcesses(this->Controller, bounds);
      this->BoundsTranslator->InitWholeBounds();
      int whole_extent[6] = {0, 8191, 0, 8191, 0, 8191};
      this->BoundsTranslator->SetWholeExtent(whole_extent);

      //      this->CachedNeuronMesh = UnstructuredGridToPolyData(vtkUnstructuredGrid::SafeDownCast(this->MeshPartitionFilter->GetOutput()), this->CachedNeuronMesh);
      this->CachedNeuronMesh = vtkPolyData::SafeDownCast(this->MeshPartitionFilter->GetOutput());
      this->MeshPartitionFilter->SetInputData((vtkPolyData*)NULL);
      //      this->MeshPartitionFilter = NULL;
      //
      redist_timer->StopTimer();
      if (this->UpdatePiece==0) {
        vtkDebugMacro(<< "ParallelRedistribution : " << redist_timer->GetElapsedTime() << " seconds");
      }
    }
#endif
    this->MeshGeneratedTime.Modified();
  }
  if (NeedToRegenerateMesh || NeedToRegernerateTime) {
    bool do_rep = (this->NumberOfTimeSteps>0) && this->GetPointArrayStatus(BBP_ARRAY_NAME_VOLTAGE);
    if (do_rep /*&& this->UpdateNumPieces==1*/) {
      try {
        this->CreateReportScalars(request, inputVector, outputVector);
      }
      catch (...) {
        vtkErrorMacro(<<"Exception caught during creation of report scalars");
      }
    }
    else {
      this->CachedNeuronMesh->GetPointData()->RemoveArray(BBP_ARRAY_NAME_VOLTAGE);
    }
    this->TimeModifiedTime.Modified();
  }

  //
  // copy internal mesh to output
  //
    if (this->ExportNeuronMesh) {
      output0->ShallowCopy(this->CachedNeuronMesh);
    }
    else if (this->ExportMorphologySkeleton) {
      output0->ShallowCopy(this->CachedMorphologySkeleton);
    }


  //
  //  output1->ShallowCopy(this->CachedMorphologySkeleton);
  //  
  load_timer->StopTimer();
  if (this->UpdatePiece==0) {
    vtkDebugMacro(<< "Mesh Load and Redistribution : " << load_timer->GetElapsedTime() << " seconds");
  }
  if (this->DeleteExperiment) {
    this->PrimaryTarget      = bbp::Target("dumy",bbp::TARGET_CELL);
    this->Partitioned_target = bbp::Target("dumy",bbp::TARGET_CELL);
    this->Experiment.clear();
    this->Microcircuit->close();
    //this->ReportReader->clearCache();
    //this->ReportMapping      
    this->OffsetMapping.clear();
  }
  return 1;
}
//-----------------------------------------------------------------------------
void vtkCircuitReader::AddOneNeuronToMesh(bbp::Neuron *neuron, const bbp::Mesh *mesh, vtkIdType Ncount, vtkPoints *points, vtkIdType *cells, vtkFieldData *field, vtkIdType &offsetN, vtkIdType &offsetC)
{  
  bool do_nrm = 1==this->GetPointArrayStatus(BBP_ARRAY_NAME_NORMAL);
  bool do_sid = 1==this->GetPointArrayStatus(BBP_ARRAY_NAME_SECTION_ID);
  bool do_nid = 1==this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONGID);
  bool do_nix = 1==this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONINDEX);
  bool do_rtn = 1==this->GetPointArrayStatus(BBP_ARRAY_NAME_RTNEURON_OPACITY);

  //
  vtkSmartPointer<vtkIntArray>   neuronId = do_nid ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_NEURONGID)) : NULL;
  vtkSmartPointer<vtkIntArray>   neuronIx = do_nix ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_NEURONINDEX)) : NULL;
  vtkSmartPointer<vtkIntArray>  sectionId = do_sid ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_SECTION_ID)) : NULL;
  vtkSmartPointer<vtkFloatArray> nvectors = do_nrm ? vtkFloatArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_NORMAL)) : NULL;
  vtkSmartPointer<vtkFloatArray> rtneuron = do_rtn ? vtkFloatArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_RTNEURON_OPACITY)) : NULL;
  //
  uint32_t         vertexCount = mesh->vertex_count();
  uint32_t         faceCount = mesh->triangle_count();
  vtkDebugMacro(<<"Neuron " << neuron->gid() << " : Triangles " << faceCount << " : Vertices " << vertexCount);
  const Vector3fs       &vertices = mesh->vertices();
  const uint16_t     *section_ids = mesh->vertex_sections().data();
  const float          *positions = mesh->vertex_relative_distances().data();
  const Vector3fs &vertex_normals = mesh->normals();
  const bbp::Morphology &morph = neuron->morphology();
  //
  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->PostMultiply();
  transform->RotateY(neuron->orientation().w()); // x,y,z,w=angle
  transform->Translate(neuron->position());      // 
  transform->Update();
  vtkIdType insertN = offsetN;
  for (uint32_t v=0 ; v<vertexCount; ++v) {
    float newPoint[3];
    transform->TransformPoint(vertices[v],newPoint); 
    points->SetPoint(offsetN, newPoint);
    uint16_t sectionID = section_ids[v];

    if (do_nrm) {
      transform->TransformNormal(vertex_normals[v],newPoint); 
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
      if (section.type() == bbp::SECTION_SOMA) {
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
      } else if (section.parent().type() == bbp::SECTION_SOMA) {
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
  const uint32_t* faces2 = mesh->triangles().data();
  for (uint32_t t=0; t<faceCount; ++t)
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
  bool do_nid = BOOL(this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONGID)); 
  bool do_nix = BOOL(this->GetPointArrayStatus(BBP_ARRAY_NAME_NEURONINDEX));
  bool do_sid = BOOL(this->GetPointArrayStatus(BBP_ARRAY_NAME_SECTION_ID));
  bool do_sty = BOOL(this->GetPointArrayStatus(BBP_ARRAY_NAME_SECTION_TYPE));
  bool do_ddr = BOOL(this->GetPointArrayStatus(BBP_ARRAY_NAME_DENDRITE_RADIUS));
  //
  vtkSmartPointer<vtkIntArray>                neuronId = do_nid ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_NEURONGID)) : NULL;
  vtkSmartPointer<vtkIntArray>                neuronIx = do_nix ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_NEURONINDEX)) : NULL;
  vtkSmartPointer<vtkIntArray>               sectionId = do_sid ? vtkIntArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_SECTION_ID)) : NULL;
  vtkSmartPointer<vtkUnsignedCharArray>    sectionType = do_sty ? vtkUnsignedCharArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_SECTION_TYPE)) : NULL;
  vtkSmartPointer<vtkFloatArray>        dendriteRadius = do_ddr ? vtkFloatArray::SafeDownCast(field->GetArray(BBP_ARRAY_NAME_DENDRITE_RADIUS)) : NULL;
  //
  const bbp::Morphology  &morph = neuron->morphology();
  const bbp::Sections &sections = morph.sections();
  const Matrix4f     &transform = neuron->global_transform();

  vtkIdType insertN = offsetN;
  vtkIdType i = 0;
  //
//  const bbp::Sections &sections = morph.neurites();
  Segments::const_iterator current_segment;
  for (Sections::const_iterator section=sections.begin(); section!=sections.end(); ++section) {
    Segments segments = section->segments();
    // allowed types are  SECTION_SOMA = 1,           //!< neuron cell body
    //                    SECTION_AXON,
    //                    SECTION_DENDRITE,           //!< general or basal dendrite (near to soma)
    //                    SECTION_APICAL_DENDRITE,    //!< apical dendrite (far from soma)
    //                    SECTION_UNDEFINED
    bbp::SectionType section_type = (*section).type();

    // if the segment has more than zero pieces, add the start point
    if (segments.begin()!=segments.end()) {
      Vector3f newPoint = transform*segments.begin()->begin().center();
      double radius = segments.begin()->begin().diameter()/2.0; 
      //
      points->SetPoint(offsetN + i, newPoint);
      if (do_sid) sectionId->SetValue(offsetN + i, section->id());
      if (do_sty) sectionType->SetValue(offsetN + i, section_type);
      if (do_ddr) dendriteRadius->SetValue(offsetN + i, radius);
      if (do_nid) neuronId->SetValue(offsetN + i, neuron->gid());
      if (do_nix) neuronIx->SetValue(offsetN + i, Ncount);
      i++;
    }
    // add one point for each segment piece
    for (Segments::const_iterator segment=segments.begin(); segment!=segments.end(); ++segment) {
      Vector3f newPoint = transform*segment->center();
      double radius = segment->diameter()/2.0; 
      //
      points->SetPoint(offsetN + i, newPoint);
      if (do_sid) sectionId->SetValue(offsetN + i, section->id());
      if (do_sty) sectionType->SetValue(offsetN + i, section_type);
      if (do_ddr) dendriteRadius->SetValue(offsetN + i, radius);
      if (do_nid) neuronId->SetValue(offsetN + i, neuron->gid());
      if (do_nix) neuronIx->SetValue(offsetN + i,Ncount);
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
    true,   // vertices
    true,   // triangles
    true,   // mapping
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
  vtkDebugMacro(<<"Allocating space for " << maxCells << " triangles");
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
  vtkDebugMacro(<<"Triangles read " << offsetC);
  vtkIdType GlobalTotalTriangles = offsetC;
  this->Controller->AllReduce(&offsetC, &GlobalTotalTriangles/*(vtkIdType*)MPI_IN_PLACE*/, 1, vtkCommunicator::SUM_OP);
  vtkDebugMacro(<<"Global number of triangles read " << GlobalTotalTriangles);

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
  //  vtkDebugMacro(<< "Report Reader Cell Target \n" << this->ReportReader->getCellTarget());

  this->ReportReader->updateMapping(this->Partitioned_target);

  // we need a cell Target object (another neuron list), so get one from the neuron list (waste of memory?)
  bbp::Cell_Target ctarget = this->Partitioned_target.cell_target();
  if (ctarget.size() != 0) {
    size_t index = 0;
    for (GIDSetCIter i=ctarget.begin(); i!=ctarget.end(); ++i, ++index) {
      this->OffsetMapping[*i] = index;
    }
  } else {
    bbp::Neurons &neurons = this->Microcircuit->neurons(); 
    size_t j = 0;
    for (bbp::Neurons::const_iterator i = neurons.begin(); i != neurons.end(); ++i,++j) {
      std::cerr << "We must provide some kind of mapping for neuron list " << std::endl;
      //            (*i)->setSimulationBufferIndex(j);
    }
  }

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
    vtkDebugMacro(<<"Allocating voltage array for N=" << this->CachedNeuronMesh->GetPoints()->GetNumberOfPoints());
    vtkIdType maxPointsN = this->NumberOfPointsBeforePartitioning;
    voltageN->Resize(maxPointsN);
    voltageN->SetNumberOfTuples(maxPointsN);
  }

  bbp::Neurons &neurons = this->Microcircuit->neurons(); 
  vtkIdType offsetN = 0, offsetM = 0;
  for (bbp::Target::cell_iterator gid=this->Partitioned_target.cell_begin(); gid!=this->Partitioned_target.cell_end(); ++gid) {
    Neurons::iterator neuron = neurons.find( *gid );
    //
    if (this->ExportMorphologySkeleton) {
      offsetM = this->AddReportScalarsToMorphologySkeleton(&*neuron, voltageM, offsetM);  
    }
    if (this->ExportNeuronMesh) {
      offsetN = this->AddReportScalarsToNeuronMesh(&*neuron, voltageN, offsetN);  
    }
  }

  if (this->ExportNeuronMesh && this->UpdateNumPieces>1 && this->ParallelRedistribution) {
    // create tmep data structures for partitioning field data
    vtkSmartPointer<vtkPointData>         input_scalars = vtkSmartPointer<vtkPointData>::New();
    vtkSmartPointer<vtkPointData>   partitioned_scalars = vtkSmartPointer<vtkPointData>::New();
    vtkSmartPointer<vtkFloatArray> partitioned_voltages = vtkSmartPointer<vtkFloatArray>::New();
    input_scalars->AddArray(voltageN);
    partitioned_voltages->SetName(voltageN->GetName());
    partitioned_scalars->AddArray(partitioned_voltages);
    partitioned_voltages->SetNumberOfTuples(this->CachedNeuronMesh->GetNumberOfPoints());
#ifdef PV_BBP_USE_ZOLTAN
    // perform the partitioning using stored partition information
    this->MeshPartitionFilter->MigratePointData(input_scalars, partitioned_scalars);
    this->CachedNeuronMesh->GetPointData()->AddArray(partitioned_voltages);
#endif
  }
  else if (this->ExportNeuronMesh) {
    this->CachedNeuronMesh->GetPointData()->AddArray(voltageN);
  }
}
//-----------------------------------------------------------------------------
vtkIdType vtkCircuitReader::AddReportScalarsToMorphologySkeleton(bbp::Neuron *neuron, vtkFloatArray *voltages, vtkIdType offsetN)
{
  const bbp::Morphology  &morph = neuron->morphology();
  const size_t index = this->OffsetMapping[neuron->gid()];
  vtkDebugMacro(<< "Neuron with GID " << neuron->gid() << " has mapping offset " << index);


  // These are the buffer offsets for this neuron 
  const brion::SectionOffsets &alloffsets = this->ReportReader->getOffsets(); // bbp:sectionoffsets
  const uint64_ts &offsets = alloffsets[index];
  
  const bbp::floatsPtr fbuffer = this->_currentFrame.getData< bbp::floatsPtr >();
  const float *buffer = fbuffer->data();
  float somaVoltage = buffer[offsets[0]];
  float rvoltage = this->RestingPotentialVoltage;

  // Finding the voltage value of the last compartment of the last axon section. 
  // This is not completely accurate but it's the best we can do to assign some color 
  // to axon sections with undefined simulation data.
  size_t cell_index = neuron->index();
  const Sections &axon = morph.axon();
  const uint16_t lastAxon = (axon.size() == 0) || 
    this->ReportReader->getCompartmentCounts()[cell_index][axon.begin()->id()] == 0 ? axon.begin()->id() : morph.soma().id();
  //
  float undefinedAxonVoltage = buffer[offsets[lastAxon] + this->ReportReader->getCompartmentCounts()[cell_index][lastAxon]];

  vtkIdType i = 0;
  const bbp::Sections &sections = morph.sections();
  Segments::const_iterator current_segment;
  for (Sections::const_iterator section=sections.begin(); section!=sections.end(); ++section) {
    Segments  segments = section->segments();
    uint16_t sectionId = section->id();

    // allowed types are  SECTION_SOMA = 1,           //!< neuron cell body
    //                    SECTION_AXON,
    //                    SECTION_DENDRITE,           //!< general or basal dendrite (near to soma)
    //                    SECTION_APICAL_DENDRITE,    //!< apical dendrite (far from soma)
    //                    SECTION_UNDEFINED
    bbp::SectionType section_type = (*section).type();

    // if the segment has more than zero pieces, add the start point
    if (segments.begin()!=segments.end()) {
      double      radius = segments.begin()->begin().diameter()/2.0; 
      //
      uint16_t compartments = this->ReportReader->getCompartmentCounts()[index][sectionId];

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
      uint16_t compartments = this->ReportReader->getCompartmentCounts()[index][sectionId];

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
return 0;
}
//-----------------------------------------------------------------------------
vtkIdType vtkCircuitReader::AddReportScalarsToNeuronMesh(bbp::Neuron *neuron, vtkFloatArray *voltages, vtkIdType offsetN)
{
  const bbp::Morphology *sdk_morph    = &neuron->morphology();
  const bbp::Mesh       *sdk_mesh     = &sdk_morph->mesh();
  //
  const uint16_ts       &section_ids = sdk_mesh->vertex_sections();
  const floats    &section_distances = sdk_mesh->vertex_relative_distances();
  //
  vtkIdType vertexCount = sdk_mesh->vertex_count();

  // get the index 
  size_t cell_index = neuron->index();
  if (cell_index==std::numeric_limits<size_t>::max()) {
    for (uint32_t i=0; i<vertexCount; ++i) {
      voltages->SetValue(offsetN + i, this->RestingPotentialVoltage);
    }
    return offsetN + vertexCount;
  }

  const size_t index = this->OffsetMapping[neuron->gid()];

  // These are the buffer offsets for this neuron 
  const brion::SectionOffsets &alloffsets = this->ReportReader->getOffsets(); // bbp:sectionoffsets
  const uint64_ts &offsets = alloffsets[index];
  
  const bbp::floatsPtr fbuffer = this->_currentFrame.getData< bbp::floatsPtr >();
  const float *buffer = fbuffer->data();
  float somaVoltage = buffer[offsets[0]];
  float rvoltage = this->RestingPotentialVoltage;

  // Finding the voltage value of the last compartment of the last axon section. 
  // This is not completely accurate but it's the best we can do to assign some color 
  // to axon sections with undefined simulation data.
  const Sections &axon = sdk_morph->axon();
  const uint16_t lastAxon = (axon.size() == 0) || 
    this->ReportReader->getCompartmentCounts()[cell_index][axon.begin()->id()] == 0 ? axon.begin()->id() : sdk_morph->soma().id();
  //
  float undefinedAxonVoltage = buffer[offsets[lastAxon] + this->ReportReader->getCompartmentCounts()[cell_index][lastAxon]];

  for (uint32_t i=0; i<vertexCount; ++i) {
    uint16_t     sectionId = section_ids[i];
    size_t num_sections = sdk_morph->size();
    uint16_t compartments = 0;
    if (sectionId < num_sections) {
      compartments = this->ReportReader->getCompartmentCounts()[cell_index][sectionId];
    }
    float position = section_distances[i];
    //
    if (compartments) {
      // Computing the relative length within section of the capsule midpoint
      uint16_t compartment = std::min(compartments - 1, (int)floor(compartments * position));
      unsigned long offset = offsets[sectionId];
      rvoltage = buffer[offset + compartment];
    } else {
      rvoltage = this->RestingPotentialVoltage;
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
  this->SILUpdateStamp++;
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
    // Get user targets for the microcircuit.
    const bbp::Targets &user_targets = this->Experiment.user_targets();
    for (bbp::Targets::const_iterator ti=user_targets.begin(); ti!=user_targets.end(); ++ti) {
      std::string name = (*ti).name();
      vtkIdType childBlock = this->SIL->AddChild(targetsRoot, childEdge);
      names.push_back(name.c_str());
      this->TargetsSelection->AddArray(name.c_str());
    }

    const bbp::Targets &default_targets = this->Experiment.targets();
    for (bbp::Targets::const_iterator ti=default_targets.begin(); ti!=default_targets.end(); ++ti) {
      std::string name = (*ti).name();
      if (!this->TargetsSelection->ArrayExists(name.c_str())) {
        vtkIdType childBlock = this->SIL->AddChild(targetsRoot, childEdge);
        names.push_back(name.c_str());
        this->TargetsSelection->AddArray(name.c_str());
      }
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
  if (this->TargetsSelection->GetArraySetting(name)!=status) {
    this->TargetsModifiedTime.Modified();
    this->Modified();
    if (status) {
      this->TargetsSelection->EnableArray(name);
    }
    else {
      this->TargetsSelection->DisableArray(name);
    }
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
