/*=========================================================================

  Project                 : Icarus
  Module                  : vtkNeuronSpikeTableSource.h

  Authors:
     John Biddiscombe     Jerome Soumagne
     biddisco@cscs.ch     soumagne@cscs.ch

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

#ifndef __vtkNeuronSpikeTableSource_h
#define __vtkNeuronSpikeTableSource_h

#include "vtkTableAlgorithm.h"
#include "vtkSmartPointer.h" // for memory safety
#include "vtkStringList.h"
// for (custom) client server arg decode define this
#define VTK_WRAPPING_CXX
#include "vtkClientServerStream.h"

#include <string>
#include <vector>
#include <unordered_map>

class vtkPoints;
class vtkCellArray;
class vtkDoubleArray;
class vtkPointData;
class vtkFloatArray;

class VTK_EXPORT vtkNeuronSpikeTableSource : public vtkTableAlgorithm
{
public:
  static vtkNeuronSpikeTableSource *New();
  vtkTypeMacro(vtkNeuronSpikeTableSource,vtkTableAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);   

  virtual void SetNameStrings(vtkStringList*);
  vtkGetObjectMacro(NameStrings, vtkStringList);

  // Description:
  // Flush will wipe any existing data so that traces can be restarted from
  // whatever time step is next supplied.
  void Flush();

  void SetDataModified();

  //
  void ClearSpikeData();
  void SetSpikeData(vtkIdType N, signed char Ids[]);
  void SetSpikeData(vtkIdType N, vtkClientServerStreamDataArg<signed char> &temp0);

protected:
   vtkNeuronSpikeTableSource();
  ~vtkNeuronSpikeTableSource();
  //
  int ProcessSpikeData(vtkTable *outdata);
  //
  int RequestInformation(
    vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector);
  //
  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);
  //
  vtkStringList  *NameStrings;
  vtkTimeStamp    DataModifiedTime;
  vtkTimeStamp    DataOpenedTime;

  // internal data variables
  int           NumberOfTimeSteps;
  //

  typedef std::tuple<uint32_t, float> spike_info_type;
  typedef std::vector<spike_info_type> spike_list_type;
  spike_list_type spikelist;
  vtkTimeStamp SpikeListModifiedTime;

  double  BinResolution;
  int     MaxTimeBins;
  double  PreviousBinStartTime;
  double  LastMaxTime;
  double  LastMinTime;

//BTX
  vtkSmartPointer<vtkPointData>   Values;
  vtkSmartPointer<vtkDoubleArray> TimeData;
  vtkSmartPointer<vtkDoubleArray> SpikeData;
  std::vector< vtkSmartPointer<vtkDoubleArray> > DataArrays;
//ETX

private:
  vtkNeuronSpikeTableSource(const vtkNeuronSpikeTableSource&);  // Not implemented.
  void operator=(const vtkNeuronSpikeTableSource&);  // Not implemented.
};

#endif
