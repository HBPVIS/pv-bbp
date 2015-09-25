/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCircuitReaderBase.cxx,v $

  Copyright (c) Sebastien Lasserre, Blue Brain Project
  All rights reserved.

=========================================================================*/

#ifndef __vtkCircuitReaderBase_h
#define __vtkCircuitReaderBase_h 
 
#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"
#include "vtkUnsignedIntArray.h"

// for (custom) client server arg decode define this
#define VTK_WRAPPING_CXX
#include "vtkClientServerStream.h"

// std::
#include <string>
#include <vector>

class vtkDataArraySelection;
class vtkMultiProcessController;
class vtkFloatArray;
class vtkMutableDirectedGraph;
class vtkMeshPartitionFilter;
class vtkParticlePartitionFilter;
class vtkUnstructuredGrid;
class vtkBoundsExtentTranslator;

// BBP-SDK
#include "BBP/Targets/Target.h"
#include "BBP/Targets/Targets.h"
#include "BBP/Datasets/compartmentReportFrame.h"
#include "BBP/Readers/compartmentReportReader.h"
#include "BBP/Experiment.h"

//----------------------------------------------------------------------------
// custom set/get macro to replace the vtk ones
//----------------------------------------------------------------------------
#undef vtkGetMacro
#undef vtkSetMacro
//
#define vtkGetMacro(name,type) \
virtual type Get##name () { \
  return this->name; \
  }

#define vtkSetMacro(name,type) \
virtual void Set##name (type _arg) \
  { \
  if (this->name != _arg) \
    { \
    this->name = _arg; \
    this->Modified(); \
    } \
  }
//
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

class VTK_EXPORT vtkCircuitReaderBase : public vtkPolyDataAlgorithm 
{
public:
  static vtkCircuitReaderBase *New();
  //vtkTypeMacro(vtkCircuitReaderBase,vtkPolyDataAlgorithm);
  vtkTypeMacro(vtkCircuitReaderBase, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of the morphology file.
  void SetFileName(char *filename);
  vtkGetStringMacro(FileName);
  
  // Description:
  // Specify the default target to load
  void SetDefaultTarget(char *target);
  vtkGetStringMacro(DefaultTarget);

  // Description:
  // Specify the report to load
  vtkSetStringMacro(ReportName);
  vtkGetStringMacro(ReportName);

  vtkSetMacro(ParallelRedistribution,int);
  vtkGetMacro(ParallelRedistribution,int);
  vtkBooleanMacro(ParallelRedistribution,int);

  vtkSetMacro(DeleteExperiment,int);
  vtkGetMacro(DeleteExperiment,int);
  vtkBooleanMacro(DeleteExperiment,int);

  // Description:
  // Set/Get the timestep to be read
  vtkSetMacro(TimeStep,int);
  vtkGetMacro(TimeStep,int);
  
  // Description:
  // Export time values as 0,1...N-1 regardless of real time values in file
  vtkSetMacro(IntegerTimeStepValues,int);
  vtkGetMacro(IntegerTimeStepValues,int);
  vtkBooleanMacro(IntegerTimeStepValues,int);

  // Description:
  // Get the number of timesteps in the file
  vtkGetMacro(NumberOfTimeSteps,int);

  vtkSetObjectMacro(SelectedGIds, vtkUnsignedIntArray);
  vtkGetObjectMacro(SelectedGIds, vtkUnsignedIntArray);

  void SetSelectedGIds(vtkIdType N, int Ids[]);
//BTX
  void SetSelectedGIds(vtkIdType N, vtkClientServerStreamDataArg<int> &temp0);
//ETX

  // Description:
  // Get the number of point or cell arrays available in the input.
  int GetNumberOfPointArrays();
  int GetNumberOfTargets();

  // Description:
  // Get the name of the point or cell array with the given index in
  // the input.
  const char* GetPointArrayName(int index);
  const char* GetTargetsName(int index);

  // Description:
  // Get/Set whether the point or cell array with the given name is to
  // be read.
  int GetPointArrayStatus(const char* name);
  int GetTargetsStatus(const char* name);

  void SetPointArrayStatus(const char* name, int status);
  void SetTargetsStatus(const char* name, int status);
  void DisableAllTargets();

  // Description:
  // Set/Get the controller use in parallel operations 
  // (set to the global controller by default)
  // If not using the default, this must be set before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

  void SetFileModified();

  // Description:
  // Every time the SIL is updated a this will return a different value.
  vtkGetMacro(SILUpdateStamp, int);

  vtkSmartPointer<vtkMutableDirectedGraph> GetSIL();

protected:
   vtkCircuitReaderBase();
  ~vtkCircuitReaderBase();

  int   FillOutputPortInformation( int port, vtkInformation* info );
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  int RequestTimeInformation(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector);

  // internally used to open report file
  int       OpenReportFile();
  //
  //
  //
  class TimeToleranceCheck: public std::binary_function<double, double, bool>
  {
  public:
    TimeToleranceCheck(double tol) { this->tolerance = tol; }
    double tolerance;
    //
    result_type operator()(first_argument_type a, second_argument_type b) const
    {
      bool result = (fabs(a-b)<=(this->tolerance));
      return (result_type)result;
    }
  };

  char         *FileName;
  char         *DefaultTarget;
  char         *ReportName;
  int           NumberOfTimeSteps;
  int           TimeStep;
  int           ActualTimeStep;
  double        TimeStepTolerance;
  double        CurrentTime;
  vtkTimeStamp  FileModifiedTime;
  vtkTimeStamp  FileOpenedTime;
  vtkTimeStamp  MeshParamsModifiedTime;
  vtkTimeStamp  MeshGeneratedTime;
  vtkTimeStamp  TimeModifiedTime;
  vtkTimeStamp  TargetsModifiedTime;
  vtkTimeStamp  InfoGeneratedTime;
  int           UpdatePiece;
  int           UpdateNumPieces;
  int           IntegerTimeStepValues;
  int           DeleteExperiment;
  int           IgnoreTime;

  // this is the selectiobn array we will get from zeq
  vtkUnsignedIntArray *SelectedGIds;
  //
//BTX
  std::vector<double> TimeStepValues;
  bbp::CompartmentReportFrame _currentFrame;
//ETX

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;
  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection *TargetsSelection;

  vtkMultiProcessController* Controller;

  bbp::Experiment                           Experiment;
  std::string                               TargetName;
  bbp::Target                               PrimaryTarget;
  bbp::Target                               Partitioned_target;
  int                                       PartitionExtents[6];
  float                          startTime;
  float                          stopTime;
  float                          timestep;
  bbp::Microcircuit_Ptr                     Microcircuit;
  bbp::CompartmentReportReaderPtr           ReportReader;
  std::map<uint32_t, size_t>  OffsetMapping;

  // 
  vtkSmartPointer<vtkMutableDirectedGraph> SIL;
  void BuildSIL();
  int SILUpdateStamp;
  int ParallelRedistribution;

  vtkIdType NumberOfPointsBeforePartitioning;

#ifdef PV_BBP_USE_ZOLTAN
  vtkSmartPointer<vtkMeshPartitionFilter>     MeshPartitionFilter;
  vtkSmartPointer<vtkParticlePartitionFilter> ParticlePartitionFilter;
  vtkSmartPointer<vtkBoundsExtentTranslator>  BoundsTranslator;
#endif
private:
  vtkCircuitReaderBase(const vtkCircuitReaderBase&); 
  void operator=(const vtkCircuitReaderBase&); 
};
 
#endif

