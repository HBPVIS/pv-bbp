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
#include "vtkCircuitReaderBase.h"
#include "BBP/Report_Specification.h"
#include "BBP/Containers/Reports_Specification.h"
//#include "SpikeData.h"

//----------------------------------------------------------------------------
#define JB_DEBUG__
//----------------------------------------------------------------------------
// extend debug macro to print the rank of the process
//----------------------------------------------------------------------------
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
//#define MANUAL_MESH_LOAD
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
using namespace bbp;
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkCircuitReaderBase);
//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkCircuitReaderBase, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
vtkCircuitReaderBase::vtkCircuitReaderBase() :
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
  this->ReportName                      = NULL;
  this->UpdatePiece                     = 0;
  this->UpdateNumPieces                 = 0;
  this->IntegerTimeStepValues           = 0;
  this->ParallelRedistribution          = 1;
  this->DeleteExperiment                = 1;
  this->IgnoreTime                      = 0;
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
#ifdef PV_BBP_USE_ZOLTAN
  this->ParticlePartitionFilter       = NULL;
  this->MeshPartitionFilter           = NULL;
  this->BoundsTranslator              = vtkSmartPointer<vtkBoundsExtentTranslator>::New();
#endif

  this->NumberOfPointsBeforePartitioning = 0;

  // when we receive a selection from zeq, this will be filled
  this->SelectedGIds = NULL;
}
//----------------------------------------------------------------------------
vtkCircuitReaderBase::~vtkCircuitReaderBase()
{
  this->SIL                           = NULL;
#ifdef PV_BBP_USE_ZOLTAN
  this->ParticlePartitionFilter       = NULL;
  this->MeshPartitionFilter           = NULL;
  this->BoundsTranslator              = NULL;
#endif
  //
  this->PointDataArraySelection->FastDelete();
  this->TargetsSelection->FastDelete();
  this->SetController(NULL);
  this->SetSelectedGIds(NULL);
  //
  delete []this->FileName;
  delete []this->DefaultTarget;
  delete []this->ReportName;
}
//----------------------------------------------------------------------------
int vtkCircuitReaderBase::FillOutputPortInformation( int port, vtkInformation* info )
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
//----------------------------------------------------------------------------
int vtkCircuitReaderBase::RequestTimeInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  int result = 1;
  // if there is no blue config supplied yet, exit quietly
  if (!this->FileName) {
    return 1;
  }
  // vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
  //  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
  //
  this->UpdatePiece = this->Controller->GetLocalProcessId(); // outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = this->Controller->GetNumberOfProcesses(); // outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  return result;
}
//----------------------------------------------------------------------------
int vtkCircuitReaderBase::RequestInformation(
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
  outInfo0->Set(CAN_HANDLE_PIECE_REQUEST(), 1);
//  outInfo0->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  bool NeedToReloadFile = (FileModifiedTime>FileOpenedTime);

  if (!vtksys::SystemTools::FileExists(this->FileName)) {
    vtkWarningMacro("File not found " << this->FileName);
    NeedToReloadFile = 0;
    result = 0;
  }
  bool NeedToRegenerateInfo = (TargetsModifiedTime>InfoGeneratedTime);

  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_NORMAL);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_NEURONGID);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_NEURONINDEX);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_SECTION_ID);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_SECTION_TYPE);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_DENDRITE_RADIUS);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_VOLTAGE);
  this->PointDataArraySelection->AddArray(BBP_ARRAY_NAME_RTNEURON_OPACITY);

  bool ok = true;
  if (NeedToReloadFile || NeedToRegenerateInfo) {
    std::string blueconfig = this->FileName;
    // -------------------------------------------------------------------   
    // Create BBP-SDK Experiment and Microcircuit to access to the neurons.
    // -------------------------------------------------------------------   
    try {
      this->Experiment.close();
      this->Experiment.open(blueconfig);
      this->Microcircuit = this->Experiment.microcircuit_ptr();
    }
    catch (std::exception &e) {
      this->Experiment.clear();
      vtkErrorMacro("An exception occurred opening the bbp::Experiment " << e.what());
    }
    //
    this->BuildSIL();
    outInfo0->Set(vtkDataObject::SIL(), this->GetSIL());
    this->FileOpenedTime.Modified();
  }

  ok = true;
  if (NeedToRegenerateInfo || NeedToReloadFile) {
    // default Target?
    this->TargetName = (this->DefaultTarget && strlen(this->DefaultTarget)>0) ? this->DefaultTarget : "";
    //
    this->PrimaryTarget = bbp::Target("empty",bbp::TARGET_CELL);
    //
    int c = 0;
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
                vtkDebugMacro(<< "Adding (user) target " << name << " to load list" );
                this->PrimaryTarget.insert(temp);
                c++;
              }
            }
            catch (std::exception &e) {
              bbp::Target temp = this->Experiment.targets().get_target(name);
              if (temp.size()>0) {
                vtkDebugMacro(<< "Adding (system) target " << name << " to list : exception " << e.what());
                this->PrimaryTarget.insert(temp);
                c++;
              }
            }
          }
          catch (std::exception &e) {
            vtkErrorMacro("Could not add target " << name << " to list : exception " << e.what());
          }
        }
      }
      // if no targets set in TargetsStatus, then use default target
      if (c==0 && this->TargetName.size()>0) {
        bbp::Target temp = this->Experiment.user_targets().get_target(this->TargetName);
        if (temp.size()>0) {
          vtkDebugMacro(<< "Adding (default) target " << this->TargetName.c_str() << " to list ");
          this->PrimaryTarget.insert(temp);
          c++;
        }
      }
    }
    catch (std::exception &e) {
      vtkErrorMacro("Could not set the target : exception " << e.what());
      this->PrimaryTarget = bbp::Target("exception",bbp::TARGET_CELL);; 
      ok = false;
    }
    this->InfoGeneratedTime.Modified();
  }

  bool needToRegenerateTimeInfo = (this->IgnoreTime==0);
  if (needToRegenerateTimeInfo) {
    // time steps of reports are in the report file
    if (ok /*&& this->UpdateNumPieces==1*/) {
      try {
        if (this->OpenReportFile()) {
          this->NumberOfTimeSteps = (this->stopTime-this->startTime)/this->timestep;
          vtkDebugMacro(<< "Number of time steps is " << this->NumberOfTimeSteps)
        }
        else {
          this->NumberOfTimeSteps = 0;
        }
      }
      catch (std::exception &e) {
    	  vtkErrorMacro("Exception caught during report file read " << e.what())
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
//    outInfo0->Set(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR(), this->BoundsTranslator);
  }
#endif

  //
  return result;
}
//----------------------------------------------------------------------------
int vtkCircuitReaderBase::OpenReportFile()
{
  std::string reportname = "";
  std::string ideal = this->ReportName ? this->ReportName : "";
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
//-----------------------------------------------------------------------------
void vtkCircuitReaderBase::BuildSIL()
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
  this->TargetsSelection->DisableAllArrays();
  TargetsModifiedTime.Modified();

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
void vtkCircuitReaderBase::SetFileName(char *filename)
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
void vtkCircuitReaderBase::SetDefaultTarget(char *target)
{
  if (this->DefaultTarget == NULL && target == NULL)
  {
    return;
  }
  if (this->DefaultTarget && target && (!strcmp(this->DefaultTarget,target)))
  {
    return;
  }
  delete [] this->DefaultTarget;
  this->DefaultTarget = NULL;

  if (target)
  {
    this->DefaultTarget = vtksys::SystemTools::DuplicateString(target);
    this->TargetsModifiedTime.Modified();
  }
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkCircuitReaderBase::SetFileModified()
{
  this->FileModifiedTime.Modified();
  this->Modified();
}

//-----------------------------------------------------------------------------
vtkSmartPointer<vtkMutableDirectedGraph> vtkCircuitReaderBase::GetSIL()
{
  return this->SIL;
}

//----------------------------------------------------------------------------
int vtkCircuitReaderBase::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
const char* vtkCircuitReaderBase::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkCircuitReaderBase::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkCircuitReaderBase::SetPointArrayStatus(const char* name, int status)
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
int vtkCircuitReaderBase::GetNumberOfTargets()
{
  return this->TargetsSelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
const char* vtkCircuitReaderBase::GetTargetsName(int index)
{
  return this->TargetsSelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkCircuitReaderBase::GetTargetsStatus(const char* name)
{
  return this->TargetsSelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkCircuitReaderBase::SetTargetsStatus(const char* name, int status)
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
void vtkCircuitReaderBase::DisableAllTargets()
{
  this->TargetsSelection->DisableAllArrays();
  this->TargetsModifiedTime.Modified();
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkCircuitReaderBase::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
    << (this->FileName ? this->FileName : "(none)") << "\n";  
}

//----------------------------------------------------------------------------
void vtkCircuitReaderBase::ClearSelectedGIds()
{
  if (!this->SelectedGIds) this->SelectedGIds = vtkUnsignedIntArray::New();
  else {
    this->SelectedGIds->Initialize();
  }
  this->MeshParamsModifiedTime.Modified();
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkCircuitReaderBase::SetSelectedGIds(vtkIdType N, unsigned int *Ids)
{
  vtkWarningMacro("SetSelectedGIds - Type 1 " << N);
  this->SelectedGIds = vtkUnsignedIntArray::New();
  this->SelectedGIds->SetArray((unsigned int*)(Ids), N, 1);
  this->SetSelectedGIds(this->SelectedGIds);
  //
  this->MeshParamsModifiedTime.Modified();
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkCircuitReaderBase::SetSelectedGIds(vtkIdType N, vtkClientServerStreamDataArg<unsigned int> &temp0)
{
  vtkWarningMacro("SetSelectedGIds - Type 2 " << N);
  unsigned int *new_data = new unsigned int[N];
  std::copy(temp0.operator unsigned int *(), temp0.operator unsigned int *()+N, new_data);
  this->SetSelectedGIds(N, new_data);
}
