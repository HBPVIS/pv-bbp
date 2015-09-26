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
 #include "vtkParticlePartitionFilter.h"
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
#include "BBP/Readers/Microcircuit_Reader.h"
#include "BBP/Readers/Mesh_Reader.h"
#include "BBP/Readers/compartmentReportReader.h"
#include "BBP/Cell_Target.h"
#include "BBP/Targets/Targets.h"
#include "BBP/Soma.h"
#include "BBP/Mesh.h"
#include "BBP/Microcircuit.h"
#include "BBP/Morphology.h"
#include "BBP/Neuron.h"
#include "BBP/Section.h"
#include "BBP/Containers/Neurons.h"
#include "BBP/Containers/Sections.h"

// Header of this Reader
#include "vtkCircuitReaderSoma.h"
#include "BBP/Report_Specification.h"
#include "BBP/Containers/Reports_Specification.h"
//#include "SpikeData.h"

//----------------------------------------------------------------------------
// override debug macro to print rank with msg
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
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
using namespace bbp;
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkCircuitReaderSoma);
//----------------------------------------------------------------------------
vtkCircuitReaderSoma::vtkCircuitReaderSoma() : vtkCircuitReaderBase()
{
  this->CachedNeuronSoma = vtkSmartPointer<vtkPolyData>::New();
  this->IgnoreTime = 1;
}
//----------------------------------------------------------------------------
vtkCircuitReaderSoma::~vtkCircuitReaderSoma()
{
  this->CachedNeuronSoma = NULL;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int vtkCircuitReaderSoma::RequestData(
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

  // get the ouptut
  vtkPointSet *output0 = vtkPointSet::SafeDownCast(outInfo0->Get(vtkDataObject::DATA_OBJECT()));

  // Which time step has been requested
  double requestedTimeValue = outInfo0->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) 
    ? outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) : 0.0; 
  //
  this->ActualTimeStep = std::find_if(
    this->TimeStepValues.begin(), this->TimeStepValues.end(),
    std::bind2nd( vtkCircuitReaderBase::TimeToleranceCheck(
    this->IntegerTimeStepValues ? 0.5 : this->TimeStepTolerance ), requestedTimeValue ))
    - this->TimeStepValues.begin();
  //
  bool NeedToRegenerateTime = false;
  if (requestedTimeValue!=this->CurrentTime) {
    this->CurrentTime = requestedTimeValue;
    NeedToRegenerateTime = true;
  }
  //
  output0->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);

  // parallel pieces info
  this->UpdatePiece = outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  vtkSmartPointer<vtkTimerLog> load_timer = vtkSmartPointer<vtkTimerLog>::New();        
  load_timer->StartTimer();

  bool NeedToRegenerateMesh =
    (FileModifiedTime>MeshGeneratedTime) ||
    (MeshParamsModifiedTime>MeshGeneratedTime) ||
    (TargetsModifiedTime>MeshGeneratedTime);
  //
  if (NeedToRegenerateMesh && this->Microcircuit) {
    //
    this->GenerateSomaPoints(request, inputVector, outputVector);

    this->NumberOfPointsBeforePartitioning = this->CachedNeuronSoma->GetNumberOfPoints();

#ifdef PV_BBP_USE_ZOLTAN
    if (this->UpdateNumPieces>1 && this->ParallelRedistribution) {
      vtkSmartPointer<vtkTimerLog> redist_timer = vtkSmartPointer<vtkTimerLog>::New();        
      redist_timer->StartTimer();
      //
      this->ParticlePartitionFilter = vtkSmartPointer<vtkParticlePartitionFilter>::New();
      this->ParticlePartitionFilter->SetInputData(this->CachedNeuronSoma);
      // tell the partition filter we can dump the input memory when needed
      this->ParticlePartitionFilter->SetInputDisposable(1);
      // for animation over time, keep the map of point send/receive
      // @TODO : Only save this if voltage reports are being loaded
      this->ParticlePartitionFilter->SetKeepInversePointLists(0);

      // release our reference count for now
      this->CachedNeuronSoma = NULL;
      // Update the partition filter with parallel information
      vtkStreamingDemandDrivenPipeline::SafeDownCast(this->ParticlePartitionFilter->GetExecutive())
        ->SetUpdateExtent(0, this->UpdatePiece, this->UpdateNumPieces, 0);
      this->ParticlePartitionFilter->Update();

      // setup bounds for future use
      this->BoundsTranslator->SetKdTree(this->ParticlePartitionFilter->GetKdtree());

      vtkDataSet *dataset = vtkDataSet::SafeDownCast(this->ParticlePartitionFilter->GetOutput());
      double bounds[6];
      dataset->GetBounds(bounds);
      this->BoundsTranslator->ExchangeBoundsForAllProcesses(this->Controller, bounds);
      this->BoundsTranslator->InitWholeBounds();
      int whole_extent[6] = {0, 8191, 0, 8191, 0, 8191};
      this->BoundsTranslator->SetWholeExtent(whole_extent);
      this->CachedNeuronSoma = vtkPolyData::SafeDownCast(this->ParticlePartitionFilter->GetOutput());
      this->ParticlePartitionFilter->SetInputData((vtkPolyData*)NULL);
      // release anything we don't need
      this->ParticlePartitionFilter = NULL;
      //
      redist_timer->StopTimer();
      if (this->UpdatePiece==0) {
        vtkDebugMacro(<< "ParallelRedistribution : " << redist_timer->GetElapsedTime() << " seconds");
      }
    }
#endif
    this->MeshGeneratedTime.Modified();
  }
  if (NeedToRegenerateMesh || NeedToRegenerateTime) {
    this->TimeModifiedTime.Modified();
  }

  //
  // copy internal mesh to output
  //
  output0->ShallowCopy(this->CachedNeuronSoma);

  load_timer->StopTimer();
  if (this->UpdatePiece==0) {
    vtkDebugMacro(<< "Mesh Load and Redistribution : " << load_timer->GetElapsedTime() << " seconds");
  }
  if (this->DeleteExperiment) {
    this->PrimaryTarget      = bbp::Target("dumy",bbp::TARGET_CELL);
    this->Partitioned_target = bbp::Target("dumy",bbp::TARGET_CELL);
    this->Experiment.clear();
    this->Microcircuit->close();
    this->OffsetMapping.clear();
  }
  return 1;
}
//-----------------------------------------------------------------------------
void vtkCircuitReaderSoma::GenerateSomaPoints(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // VTK arrays
  vtkSmartPointer<vtkPoints>                    points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray>                  verts = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPointData>              pointdata = vtkSmartPointer<vtkPointData>::New();
  vtkSmartPointer<vtkFloatArray>        dendriteRadius;
  vtkSmartPointer<vtkUnsignedCharArray>    sectionType;
  vtkSmartPointer<vtkIntArray>               sectionId;
  vtkSmartPointer<vtkIntArray>                neuronId;
  vtkSmartPointer<vtkIntArray>                neuronIx;

  // Load just neurons for this target so we can partition them
  try {
    if (this->Microcircuit->neurons().size()==0) {
      this->Microcircuit->load(this->PrimaryTarget, bbp::NEURONS);
    }
  }
  catch (std::exception &e) {
    vtkErrorMacro("Caught an exception during Microcircuit->load " << e.what());
  }

  //
  bbp::Neurons &neurons = this->Microcircuit->neurons();
  //
  int WholeExtent[6] = { 0, static_cast<int>(neurons.size()), 0, 0, 0, 0 };
  if (this->SelectedGIds!=NULL) {
    WholeExtent[1] = (size_t)(this->SelectedGIds->GetNumberOfTuples());
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

  // create a new target for our subrange of neurons, clear any contests first.
  this->Partitioned_target = bbp::Target("ParaViewCells", bbp::TARGET_CELL);
  bbp::Target NeuronsToLoad_target = bbp::Target("Temporary", bbp::TARGET_CELL);

  // if the user has passed a GId array (selection), we should use them directly
  // onyl add neurons to our target if we have not loaded them before
  if (this->SelectedGIds!=NULL) {
    unsigned int *raw_data = this->SelectedGIds->GetPointer(0);
    for (int i=PartitionExtents[0]; i<PartitionExtents[1]; i++) {
      if (soma_map.find(raw_data[i])==soma_map.end()) {
        NeuronsToLoad_target.insert(raw_data[i]);
      }
      this->Partitioned_target.insert(raw_data[i]);
    }
  }
  else {
    Neurons::iterator neuron = neurons.begin();
    std::advance(neuron, PartitionExtents[0]);
    for (int i=PartitionExtents[0]; i<PartitionExtents[1]; ++i, ++neuron) {
      if (soma_map.find(neuron->gid())==soma_map.end()) {
        NeuronsToLoad_target.insert(neuron->gid());
      }
      this->Partitioned_target.insert(neuron->gid());
    }
  }

  //
  // reserve space for coordinates
  //
  vtkIdType maxPoints = PartitionExtents[1]-PartitionExtents[0];
  vtkIdType maxCells  = maxPoints;
  points->GetData()->Resize(maxPoints);
  points->SetNumberOfPoints(maxPoints);
  //
  // reserve space for cells
  //
  vtkIdType *cells = verts->WritePointer(maxCells, 2*(maxCells));
  //
  // reserve space for each scalar/vector fields
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

  // Load neurons for the target so we can partition them
  // we need morphologies to get radius of soma etc
  try {
    if (!NeuronsToLoad_target.empty()) {
      this->Microcircuit->load(NeuronsToLoad_target, bbp::NEURONS | bbp::MORPHOLOGIES);
    }
  }
  catch (std::exception &e) {
    vtkErrorMacro("Caught an exception during Microcircuit->load " << e.what());
  }

  //
  // loop over partitioned neurons
  //
  int      index = 0;
  double   radius;
  Vector3f newPoint;
  int      section_Id;
  unsigned char section_Type;
  //
  for (bbp::Target::cell_iterator gid=this->Partitioned_target.cell_begin(); gid!=this->Partitioned_target.cell_end(); ++gid, ++index) {
    soma_map_type::iterator it = soma_map.find(*gid);
    if (it==soma_map.end()) {
      Neurons::iterator neuron = neurons.find( *gid );
      radius          = neuron->soma().mean_radius();
      newPoint        = neuron->soma().position();
      section_Id      = neuron->sections().begin()->id();
      section_Type    = neuron->sections().begin()->type();
      soma_map[*gid]  = std::make_tuple(radius, newPoint, section_Id, section_Type);
      neuron->clear();
    }
    else {
      radius       = std::get<0>(it->second);
      newPoint     = std::get<1>(it->second);
      section_Id   = std::get<2>(it->second);
      section_Type = std::get<3>(it->second);
    }
    //
    points->SetPoint(index, newPoint);
    if (do_sid) sectionId->SetValue(index, section_Id);
    if (do_sty) sectionType->SetValue(index, section_Type);
    if (do_ddr) dendriteRadius->SetValue(index, radius);
    if (do_nid) neuronId->SetValue(index, *gid);
    if (do_nix) neuronIx->SetValue(index, index + PartitionExtents[0]);
    //
    cells[index*2]   = 1;
    cells[index*2+1] = index;
  }
  //
  this->CachedNeuronSoma->SetPoints(points);
  this->CachedNeuronSoma->SetVerts(verts);
  this->CachedNeuronSoma->GetPointData()->ShallowCopy(pointdata);
  vtkDebugMacro("Completed GenerateSomaPoints");
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
