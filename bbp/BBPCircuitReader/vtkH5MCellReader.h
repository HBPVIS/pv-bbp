/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkH5MCellReader.h
  Revision of last commit : $Rev: 754 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2009-01-09 13:40:38 +0100 #$

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing 
  1) This copyright notice appears on all copies of source code 
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it 
  must not be reformatted such that the indentation, bracketing or 
  overall style is modified significantly. 

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
// .NAME vtkH5MCellReader - Write H5Part (HDF5) Particle files
// .SECTION Description
// vtkH5MCellReader reads compatible with H5Part : documented here
// http://amas.web.psi.ch/docs/H5Part-doc/h5part.html 

#ifndef __vtkH5MCellReader_h
#define __vtkH5MCellReader_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"
#include <string>
#include <vector>
#include "hdf5.h"

class vtkDataArraySelection;
class vtkMultiProcessController;
class vtkBoundsExtentTranslator;

class VTK_EXPORT vtkH5MCellReader : public vtkPolyDataAlgorithm
{
public:
  static vtkH5MCellReader *New();
  vtkTypeMacro(vtkH5MCellReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);   

  // Description:
  // Specify file name.
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
  // When set (default no), the reader will generate a vertex cell
  // for each point/particle read. When using the points directly
  // this is unnecessary and time can be saved by omitting cell generation
  // vtkPointSpriteMapper does not require them.
  // When using ParaView, cell generation is recommended, without them
  // many filter operations are unavailable
  vtkSetMacro(GenerateVertexCells, int);
  vtkGetMacro(GenerateVertexCells, int);
  vtkBooleanMacro(GenerateVertexCells, int);

  // Description:
  // Normally, a request for data at time t=x, where x is either before the start of
  // time for the data, or after the end, will result in the first or last
  // timestep of data to be retrieved (time is clamped to max/min values).
  // Forsome applications/animations, it may be desirable to not display data
  // for invalid times. When MaskOutOfTimeRangeOutput is set to ON, the reader
  // will return an empty dataset for out of range requests. This helps
  // avoid corruption of animations.
  vtkSetMacro(MaskOutOfTimeRangeOutput, int);
  vtkGetMacro(MaskOutOfTimeRangeOutput, int);
  vtkBooleanMacro(MaskOutOfTimeRangeOutput, int);

  // Description:
  // An H5Part file may contain multiple arrays
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
  // Set/Get the controller use in compositing (set to
  // the global controller by default)
  // If not using the default, this must be called before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

  void SetFileModified();

protected:
   vtkH5MCellReader();
  ~vtkH5MCellReader();
  //
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  //
  virtual int  OpenFile();
  virtual void CloseFile();
  
  // Under normal circumstances, when loading data and animating though timesteps
  // one does not want to close the file between steps (calls to ExecuteInfo/Data)
  // but subclasses (e.g. dsm) do need to close the file for real. We therefore
  // call CloseFileIntermediate when we can leave it open and subclasses can decide
  // to act on it, or do nothing. By default, do nothing.
  virtual void CloseFileIntermediate();
  vtkSmartPointer<vtkDataArray> ReadDataSetIntoDataArray(char *name, vtkIdType s, vtkIdType e);

  //
  // Internal Variables
  //
  char         *FileName;
  vtkSmartPointer<vtkDataArray> offsets;
  int           NumberOfTimeSteps;
  int           TimeStep;
  int           ActualTimeStep;
  double        TimeStepTolerance;
  int           GenerateVertexCells;
  hid_t         H5FileId;
  vtkTimeStamp  FileModifiedTime;
  vtkTimeStamp  FileOpenedTime;
  int           UpdatePiece;
  int           UpdateNumPieces;
  int           MaskOutOfTimeRangeOutput;
  int           TimeOutOfRange;
  int           IntegerTimeStepValues;
  //
//BTX
  std::vector<double> TimeStepValues;
//ETX

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection     *PointDataArraySelection;
  vtkMultiProcessController *Controller;

private:
  vtkH5MCellReader(const vtkH5MCellReader&);  // Not implemented.
  void operator=(const vtkH5MCellReader&);  // Not implemented.
};

#endif
