/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCircuitReader.cxx,v $

  Copyright (c) Sebastien Lasserre, Blue Brain Project
  All rights reserved.

=========================================================================*/

#ifndef __vtkCircuitReader_h
#define __vtkCircuitReader_h 
 
#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"
#include "vtkUnsignedIntArray.h"
#define VTK_WRAPPING_CXX
#include "vtkClientServerStream.h"
#include <string>
#include <vector>


class vtkDataArraySelection;
class vtkMultiProcessController;
class vtkFloatArray;
class vtkMutableDirectedGraph;
class vtkMeshPartitionFilter;
class vtkUnstructuredGrid;
class vtkBoundsExtentTranslator;
// BBP-SDK

//#include "BBP/BBP.h"
#include "BBP/Targets/Target.h"
#include "BBP/Targets/Targets.h"
#include "BBP/Datasets/compartmentReportFrame.h"
#include "BBP/Readers/compartmentReportReader.h"
#include "BBP/Experiment.h"

//
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
class VTK_EXPORT vtkCircuitReader : public vtkPolyDataAlgorithm 
{
public:
  static vtkCircuitReader *New();
  //vtkTypeMacro(vtkCircuitReader,vtkPolyDataAlgorithm);
  vtkTypeMacro(vtkCircuitReader,vtkPolyDataAlgorithm);
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

  void SetExportMorphologySkeleton(int n) {
    if (this->ExportMorphologySkeleton != n) { 
      this->ExportMorphologySkeleton = n; 
      this->MeshParamsModifiedTime.Modified();
      this->Modified(); 
    } 
  }
  vtkGetMacro(ExportMorphologySkeleton,int);
  vtkBooleanMacro(ExportMorphologySkeleton,int);

  void SetExportNeuronMesh(int n) {
    if (this->ExportNeuronMesh != n) { 
      this->ExportNeuronMesh = n; 
      this->MeshParamsModifiedTime.Modified();
      this->Modified(); 
    } 
  }
  vtkGetMacro(ExportNeuronMesh,int);
  vtkBooleanMacro(ExportNeuronMesh,int);

  vtkSetMacro(ParallelRedistribution,int);
  vtkGetMacro(ParallelRedistribution,int);
  vtkBooleanMacro(ParallelRedistribution,int);

  void SetMaximumNumberOfNeurons(int n) {
    if (this->MaximumNumberOfNeurons != n) { 
      this->MaximumNumberOfNeurons = n; 
      this->MeshParamsModifiedTime.Modified();
      this->Modified(); 
    } 
  }
  vtkGetMacro(MaximumNumberOfNeurons,int);

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

  vtkSetObjectMacro(SelectedGIds, vtkUnsignedIntArray);
  vtkGetObjectMacro(SelectedGIds, vtkUnsignedIntArray);

  void SetSelectedGIds(vtkIdType N, int Ids[]);
//BTX
  void SetSelectedGIds(vtkIdType N, vtkClientServerStreamDataArg<int> &temp0);
//ETX

  // Description:
  // Get the number of timesteps in the file
  vtkGetMacro(NumberOfTimeSteps,int);

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

  // Description:
  // hyperpolarized (voltage near -85 mV)
  vtkSetMacro(HyperPolarizedVoltage,double);
  vtkGetMacro(HyperPolarizedVoltage,double);

  // Description:
  // depolarized (voltage near -50 mV)
  vtkSetMacro(DePolarizedVoltage,double);
  vtkGetMacro(DePolarizedVoltage,double);

  // Description:
  // resting potential (voltage near -65 mV)
  vtkSetMacro(RestingPotentialVoltage,double);
  vtkGetMacro(RestingPotentialVoltage,double);
  
protected:
   vtkCircuitReader();
  ~vtkCircuitReader();

  int   FillOutputPortInformation( int port, vtkInformation* info );
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
  int RequestTimeInformation(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector);

  // top level functions which generate meshes/scalars for the whole data
  void GenerateNeuronMesh(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  void GenerateMorphologySkeleton(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  void CreateReportScalars(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  // generate mesh/scalars on a per neuron basis
  void      AddOneNeuronToMesh(bbp::Neuron *neuron, const bbp::Mesh *mesh, vtkIdType Ncount, vtkPoints *points, vtkIdType *cells, vtkFieldData *field, vtkIdType &offsetN, vtkIdType &offsetC);
  void      AddOneNeuronToMorphologySkeleton(bbp::Neuron *neuron, vtkIdType Ncount, vtkPoints *points, vtkIdType *cells, vtkFieldData *field, vtkIdType &offsetN, vtkIdType &offsetC);
  vtkIdType AddReportScalarsToMorphologySkeleton(bbp::Neuron *neuron, vtkFloatArray *voltages, vtkIdType offsetN);
  vtkIdType AddReportScalarsToNeuronMesh(bbp::Neuron *neuron, vtkFloatArray *voltages, vtkIdType offsetN);

  // internally used to open report file
  int       OpenReportFile();
  //
  //
  //

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
  int ExportMorphologySkeleton;
  int ExportNeuronMesh;
  int ParallelRedistribution;
  int MaximumNumberOfNeurons;

  double HyperPolarizedVoltage;
  double DePolarizedVoltage;
  double RestingPotentialVoltage;

  vtkIdType NumberOfPointsBeforePartitioning;

  vtkSmartPointer<vtkPolyData>               CachedNeuronMesh;
  vtkSmartPointer<vtkPolyData>               CachedMorphologySkeleton;
#ifdef PV_BBP_USE_ZOLTAN
  vtkSmartPointer<vtkMeshPartitionFilter>    MeshPartitionFilter;
  vtkSmartPointer<vtkBoundsExtentTranslator> BoundsTranslator;
#endif
private:
  vtkCircuitReader(const vtkCircuitReader&); 
  void operator=(const vtkCircuitReader&); 
};
 
#endif

