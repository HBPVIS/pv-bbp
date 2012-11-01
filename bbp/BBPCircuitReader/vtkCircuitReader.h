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
#include <vtkstd/string>
#include <vtkstd/vector>

class vtkDataArraySelection;
class vtkMultiProcessController;
class vtkBoundsExtentTranslator;

// BBP-SDK

#include <BBP/common.h>
#include "BBP/Microcircuit/Experiment.h"
 
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
  // An file may contain multiple arrays
  // a GUI (eg Paraview) can provide a mechanism for selecting which data arrays
  // are to be read from the file. The PointArray variables and members can
  // be used to query the names and number of arrays available
  // and set the status (on/off) for each array, thereby controlling which
  // should be read from the file. Paraview queries these point arrays after
  // the (update) information part of the pipeline has been updated, and before the
  // (update) data part is updated.
  int         GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int         GetPointArrayStatus(const char* name);
  void        SetPointArrayStatus(const char* name, int status);
  void        DisableAll();
  void        EnableAll();
  void        Disable(const char* name);
  void        Enable(const char* name);
  //
  int         GetNumberOfPointArrayStatusArrays() { return GetNumberOfPointArrays(); }
  const char* GetPointArrayStatusArrayName(int index) { return GetPointArrayName(index); }
  int         GetPointArrayStatusArrayStatus(const char* name) { return GetPointArrayStatus(name); }
  void        SetPointArrayStatusArrayStatus(const char* name, int status) { SetPointArrayStatus(name, status); }

  // Description:
  // Set/Get the controller use in parallel operations 
  // (set to the global controller by default)
  // If not using the default, this must be set before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

  void SetFileModified();


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
  typedef vtkstd::vector<vtkstd::string>  stringlist;
  vtkstd::vector<stringlist>              FieldArrays;
  // For Bounding boxes if present
//ETX

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;

  vtkMultiProcessController* Controller;

    bbp::Experiment experiment;
 
private:
  vtkCircuitReader(const vtkCircuitReader&); 
  void operator=(const vtkCircuitReader&); 
};
 
#endif

