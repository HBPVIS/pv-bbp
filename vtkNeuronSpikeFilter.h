/*=========================================================================

Program:   Visualization Toolkit
Module:    $RCSfile: vtkNeuronSpikeFilter.h,v $

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkNeuronSpikeFilter - Create a pseudo moving average of data using an Exponential Decay Function
// .SECTION Description
// 

#ifndef __vtkNeuronSpikeFilter_h
#define __vtkNeuronSpikeFilter_h

#include "vtkExponentialDecayFilter.h"
#include "vtkSmartPointer.h" // required
#include <vector>     // required
#include <unordered_map>     // required

// for (custom) client server arg decode define this
#define VTK_WRAPPING_CXX
#include "vtkClientServerStream.h"

class vtkDataSet;
class vtkDataArraySelection;
class vtkPointData;
class vtkFloatArray;

class VTK_EXPORT vtkNeuronSpikeFilter : public vtkExponentialDecayFilter
{
public:
  static vtkNeuronSpikeFilter *New();
  vtkTypeMacro(vtkNeuronSpikeFilter, vtkExponentialDecayFilter);

  void ClearSpikeData();
  void SetSpikeData(vtkIdType N, signed char Ids[]);
  void SetSpikeData(vtkIdType N, vtkClientServerStreamDataArg<signed char> &temp0);

  vtkSetMacro(SpikeThreshold, int);
  vtkGetMacro(SpikeThreshold, int);

  int BuildGIdIndexMap(vtkPointSet *inData);
  int ProcessSpikeData(vtkPointSet *outdata);

  int RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector);
  
  typedef std::tuple<uint32_t, float>  spike_info_type;
  typedef std::vector<spike_info_type> spike_list_type;
  typedef std::vector<vtkIdType>       pointIdList;
  spike_list_type spikelist;

  typedef std::unordered_map<uint32_t, pointIdList> gid_index_map_type;
  int                 NeedToRegenerateMap;
  gid_index_map_type  gid_index_map;
  vtkTimeStamp        gid_index_map_time;
  vtkTimeStamp        SpikeListModifiedTime;
  vtkSmartPointer<vtkFloatArray> last_spikevalue;

  int UpdatePiece;
  int UpdateNumPieces;
  int SpikeThreshold;
protected:
   vtkNeuronSpikeFilter();
  ~vtkNeuronSpikeFilter();
/*
  virtual int ExecuteInformation(vtkInformation *,
    vtkInformationVector **,
    vtkInformationVector *);

  virtual int RequestData(vtkInformation *,
    vtkInformationVector **,
    vtkInformationVector *);

  // Description:
  virtual vtkDataSet *DecayDataSet(vtkDataSet *in1, double timevalue);

  // Description:
  virtual vtkDataArray *DecayDataArray(double timevalue, vtkDataArray **arrays, vtkIdType N);
*/

/*
  // internal variables
  bool   FirstIteration;
  double HighFrequencyDelta;
  int    HighFrequencyResponse;
  int    OutputAbsoluteValue;
  int    ClampAndNormalizeOutput;
  double NormalizedRange[2];
  double DecayFactor;
  char  *ArrayNamePrefix;

  double TimeStepTolerance;
  double LastUpdateTime;

  #undef VTK_WRAPPING_CXX
  //BTX
  std::vector<double>           TimeStepValues;
  vtkSmartPointer<vtkPointData> LastPointData;
  //ETX

  // To allow paraview gui to enable/disable scalars
  vtkSmartPointer<vtkDataArraySelection> PointDataArraySelection;
*/
private:
  vtkNeuronSpikeFilter(const vtkNeuronSpikeFilter&);  // Not implemented.
  void operator=(const vtkNeuronSpikeFilter&);  // Not implemented.
};

#endif
