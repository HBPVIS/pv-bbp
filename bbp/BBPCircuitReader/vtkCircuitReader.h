/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCircuitReader.cxx,v $

  Copyright (c) Sebastien Lasserre, Blue Brain Project
  All rights reserved.

=========================================================================*/

#ifndef __vtkCircuitReader_h
#define __vtkCircuitReader_h 
 
#include "vtkPolyDataAlgorithm.h"
#include "vtkBoundingBox.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkSmartPointer.h"
#include <string>
#include <vector>


class vtkDataArraySelection;
class vtkMultiProcessController;

// BBP-SDK

#include <BBP/common.h>
#include "BBP/Microcircuit/Morphology.h"
#include "BBP/Microcircuit/Experiment.h"
#include "BBP/Common/Math/Geometry/Rotation.h "

class VTK_EXPORT vtkCircuitReader : public vtkPolyDataAlgorithm 
{
public:
  static vtkCircuitReader *New();
  //vtkTypeMacro(vtkCircuitReader,vtkPolyDataAlgorithm);
  vtkTypeRevisionMacro(vtkCircuitReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  // Description:
  // Specify file name of the morphology file.
  void SetFileName(char *filename);
  vtkGetStringMacro(FileName);
  
  vtkSetMacro(MeshType,int);
  vtkGetMacro(MeshType,int);

  vtkSetMacro(GenerateNormalVectors,int);
  vtkGetMacro(GenerateNormalVectors,int);
  vtkBooleanMacro(GenerateNormalVectors,int);

  vtkSetMacro(Random,int);
  vtkGetMacro(Random,int);
  vtkBooleanMacro(Random,int);

  vtkSetMacro(MaximumNumberOfNeurons,int);
  vtkGetMacro(MaximumNumberOfNeurons,int);

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

  void CreateDatasetFromMorphology(const bbp::Morphology *morph, vtkPoints *points, vtkCellArray *lines, vtkFieldData *field, const bbp::Transform_3D<bbp::Micron> &transform);

protected:
  vtkCircuitReader();
  ~vtkCircuitReader();

  int   FillOutputPortInformation( int port, vtkInformation* info );
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  //
  char         *FileName;
  int           NumberOfTimeSteps;
  int           TimeStep;
  int           ActualTimeStep;
  double        TimeStepTolerance;
  vtkTimeStamp  FileModifiedTime;
  vtkTimeStamp  FileOpenedTime;
  int           UpdatePiece;
  int           UpdateNumPieces;
  int           IntegerTimeStepValues;
  //
//BTX
  vtkstd::vector<double>                  TimeStepValues;
  // For Bounding boxes if present
//ETX

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;
  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection *TargetsSelection;

  vtkMultiProcessController* Controller;

  bbp::Experiment experiment;

  // 
  vtkSmartPointer<vtkMutableDirectedGraph> SIL;
  void BuildSIL();
  int SILUpdateStamp;
  int MeshType;
  int GenerateNormalVectors;
  int Random;
  int MaximumNumberOfNeurons;

private:
  vtkCircuitReader(const vtkCircuitReader&); 
  void operator=(const vtkCircuitReader&); 
};
 
#endif

