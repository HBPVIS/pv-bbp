/*=========================================================================

Program:   Visualization Toolkit
Module:    $RCSfile: vtkNeuronSpikeFilter.cxx,v $

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkNeuronSpikeFilter.h"

#include "vtkCellData.h"
#include "vtkCompositeDataIterator.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkSmartPointer.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataArraySelection.h"
#include "vtkMath.h"
#include "vtkTimerLog.h"
//
#include "vtkCircuitReaderBase.h"
//
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkNeuronSpikeFilter);
//----------------------------------------------------------------------------
double decaycompute(double avg, double decay, double lasttime, double thistime, double value) 
{
  double alpha = 1.0 - 1.0/exp(decay*abs(thistime-lasttime));
  avg += alpha*(value - avg);
  return avg;
}
//----------------------------------------------------------------------------
vtkNeuronSpikeFilter::vtkNeuronSpikeFilter()
{
  this->NeedToRegenerateMap = 0;
  this->last_spikevalue = vtkSmartPointer<vtkFloatArray>::New();

/*
  this->DecayFactor               = 100.0;
  this->HighFrequencyResponse     = 0;
  this->HighFrequencyDelta        = 20.0;
  this->ArrayNamePrefix           = NULL;
  this->LastPointData             = vtkSmartPointer<vtkPointData>::New();
  this->LastUpdateTime            = 0.0;
  this->FirstIteration            = true;
  this->OutputAbsoluteValue       = 1;
  this->ClampAndNormalizeOutput   = 1;
  this->NormalizedRange[0] = -2;
  this->NormalizedRange[1] =  2;

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->PointDataArraySelection   = vtkSmartPointer<vtkDataArraySelection>::New();
*/
}
//----------------------------------------------------------------------------
vtkNeuronSpikeFilter::~vtkNeuronSpikeFilter()
{
//  this->LastPointData = NULL;
}

//----------------------------------------------------------------------------
void vtkNeuronSpikeFilter::ClearSpikeData()
{
}

//----------------------------------------------------------------------------
void vtkNeuronSpikeFilter::SetSpikeData(vtkIdType N, signed char Ids[])
{
  vtkWarningMacro("SetSpikeData - Type 1 " << N);
}

//----------------------------------------------------------------------------
void vtkNeuronSpikeFilter::SetSpikeData(vtkIdType N, vtkClientServerStreamDataArg<signed char> &temp0)
{
  vtkWarningMacro("SetSpikeData - Type 2 " << N);
  //
  spike_info_type *begin = (spike_info_type *) (temp0.operator signed char *());
  this->spikelist.assign(begin, begin + N);

  this->Modified();
  this->SpikeListModifiedTime.Modified();
}
//----------------------------------------------------------------------------
int vtkNeuronSpikeFilter::BuildGIdIndeMap(vtkPointSet *inData)
{
  vtkPointData *pd = inData->GetPointData();
  if (!pd) {
    vtkErrorMacro("Cannot process spikes without GIds");
    return 0;
  }
  vtkDataArray *gids = pd->GetArray(BBP_ARRAY_NAME_NEURONGID);
  vtkIntArray *neuronGId = vtkIntArray::SafeDownCast(gids);
  if (!gids || !neuronGId) {
    vtkErrorMacro("Cannot process spikes without GIds");
    return 0;
  }
  vtkIdType N = neuronGId->GetNumberOfTuples();
  //
  this->gid_index_map.clear();
  this->gid_index_map.reserve(N);
  //
  for (vtkIdType i=0; i<N; ++i) {
    uint32_t Gid = neuronGId->GetValue(i);
    this->gid_index_map[Gid] = i;
  }
  this->gid_index_map_time.Modified();
  return 1;
}

//----------------------------------------------------------------------------
int vtkNeuronSpikeFilter::ProcessSpikeData(vtkPointSet *outdata)
{
  vtkPointData *pd = outdata->GetPointData();
  vtkIdType N = pd->GetNumberOfTuples();
  //
  vtkSmartPointer<vtkFloatArray> spikevalue = vtkSmartPointer<vtkFloatArray>::New();
  spikevalue->SetName("SpikeValue");
  spikevalue->SetNumberOfTuples(N);
  float *spike_data = spikevalue->GetPointer(0);
  float *last_spike_data = last_spikevalue->GetPointer(0);
  if (last_spikevalue->GetNumberOfTuples()!=N) {
    last_spikevalue->SetNumberOfTuples(N);
    last_spike_data = last_spikevalue->GetPointer(0);
    // clear data
    for (vtkIdType i=0; i<N; ++i) {
      last_spike_data[i] = 0;
    }
  }

  // clear data
  for (vtkIdType i=0; i<pd->GetNumberOfTuples(); ++i) {
    spike_data[i] = 0;
  }
  //
  for (const spike_info_type &it : spikelist) {
    gid_index_map_type::iterator g_it = gid_index_map.find(std::get<0>(it));
    if (g_it != gid_index_map.end()) {
      vtkIdType index = std::get<1>(*g_it);
      spike_data[index] = last_spike_data[index] + 1.0;
    }
  }
  //
  last_spikevalue = spikevalue;
  pd->AddArray(spikevalue);

  return 1;
}

//----------------------------------------------------------------------------
int vtkNeuronSpikeFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo0 = outputVector->GetInformationObject(0);

  vtkPointSet *output0 = vtkPointSet::SafeDownCast(outInfo0->Get(vtkDataObject::DATA_OBJECT()));

  vtkPointSet *data0 = vtkPointSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Which time step has been requested
  double requestedTimeValue = outInfo0->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) 
    ? outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) : 0.0; 

  //
  output0->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);

  // parallel pieces info
  this->UpdatePiece = outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  vtkSmartPointer<vtkTimerLog> load_timer = vtkSmartPointer<vtkTimerLog>::New();        
  load_timer->StartTimer();

  bool NeedToRegenerateMap =
    (data0->GetMTime( ) > gid_index_map_time);
  //
  if (NeedToRegenerateMap) {
    this->BuildGIdIndeMap(data0);
  }

  output0->ShallowCopy(data0);
  this->ProcessSpikeData(output0);

//  vtkPointData *pd = data0->GetPointData();
//  vtkSmartPointer<vtkPointData> new_pd = vtkSmartPointer<vtkPointData>::New();
//  new_pd->ShallowCopy(pd);
  // copy the onput to the output


//  if (data0->GetInformation()->Has(vtkDataObject::DATA_GEOMETRY_UNMODIFIED()))
//  {
//    outData->GetInformation()->Set(vtkDataObject::DATA_GEOMETRY_UNMODIFIED(),1);
//  }

  outInfo0->Set(vtkDataObject::DATA_OBJECT(),output0);

  // stamp this new dataset with a time key
  output0->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),requestedTimeValue);
  //
//  this->LastUpdateTime = upTime;
//  this->LastPointData->ShallowCopy(outData->GetPointData());
  return 1;
}

/*
//----------------------------------------------------------------------------
// Change the information
int vtkNeuronSpikeFilter::ExecuteInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  vtkDataSet *inData = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPointData *pd = inData->GetPointData();

  //
  // when information is regenerated, the arrays might change (bah!)
  //
  if (pd->GetNumberOfArrays()>0) {
    vtkSmartPointer<vtkDataArraySelection> TempArraySelection = vtkSmartPointer<vtkDataArraySelection>::New();
    for (int i=0; i<pd->GetNumberOfArrays(); i++) {
      const char *name = pd->GetArray(i)->GetName();
      TempArraySelection->AddArray(name);
      if (this->PointDataArraySelection->ArrayIsEnabled(name)) {
        TempArraySelection->EnableArray(name);
      }
      else {
        TempArraySelection->DisableArray(name);
      }
    }
    this->PointDataArraySelection->CopySelections(TempArraySelection);
  }

  double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  if (upTime<=this->LastUpdateTime && !this->FirstIteration) {
    this->FirstIteration = true;
    this->LastPointData = vtkSmartPointer<vtkPointData>::New();
  }
  return 1;
}
//----------------------------------------------------------------------------
// templated difference function
template <class T>
void vtkExponentialDecayCompute(vtkNeuronSpikeFilter *tdf,
  bool first, double decay, double lasttime, double thistime,
                               vtkDataArray *output,
                               vtkDataArray **arrays,
                               vtkIdType numComp,                                    
                               T *)
{
  T *outData = static_cast<T*>(output->GetVoidPointer(0));
  T *inData0 = static_cast<T*>(arrays[0]->GetVoidPointer(0));
  T *inData1 = NULL;
  if (arrays[1]) {
    inData1 = static_cast<T*>(arrays[1]->GetVoidPointer(0));
  }
  bool    HFR = tdf->GetHighFrequencyResponse();
  double  HFD = tdf->GetHighFrequencyDelta();
  bool    ABS = tdf->GetOutputAbsoluteValue();
  bool    NOR = tdf->GetClampAndNormalizeOutput();
  double *NRM = tdf->GetNormalizedRange();
  //
  vtkIdType N = arrays[0]->GetNumberOfTuples();
  for (vtkIdType t=0; t<N; ++t)
  {
    T *value = &inData0[t*numComp];
    double vv = *value;
    double pv = 0.0;
    if (arrays[1]) {
      inData1 = static_cast<T*>(arrays[1]->GetVoidPointer(0));
      pv = inData1[t*numComp];
    }
    for (int c=0; c<numComp; ++c) {
      double temp;
      if (first || (HFR && abs(vv)>HFD)) {
        temp = ABS ? (vv) : static_cast<T>(vv);
      }
      else {
        temp = decaycompute(pv, decay, lasttime, thistime, vv);
        temp = ABS ? abs(temp) : temp;
      }
      if (NRM) {
        *outData++ = vtkMath::ClampAndNormalizeValue(temp, NRM);
      }
      else {
        *outData++ = static_cast<T>(temp);
      }
    }
  }
  output->SetNumberOfTuples(N);
}
//----------------------------------------------------------------------------
vtkDataArray *vtkEDFNewArray(vtkDataArray *da, vtkIdType Nc, vtkIdType Nt, const char *prefix)
{
  //
  // Create the array
  //
  vtkAbstractArray *aa = da->CreateArray(da->GetDataType());
  vtkDataArray *output = vtkDataArray::SafeDownCast(aa);
  //
  // initialize 
  //
  output->SetNumberOfComponents(Nc);
  output->SetNumberOfTuples(Nt);
  std::string newname = std::string(prefix) + da->GetName();
  output->SetName(newname.c_str());
  return output;
}  
//----------------------------------------------------------------------------
vtkDataArray *vtkNeuronSpikeFilter::DecayDataArray(double timevalue, vtkDataArray **arrays, vtkIdType Nt)
{
  //
  // Create the output array
  //
  int Nc = arrays[0]->GetNumberOfComponents();
  vtkDataArray *output = vtkEDFNewArray(arrays[0], Nc, Nt, this->ArrayNamePrefix);

  // now do the interpolation
  switch (arrays[0]->GetDataType())
  {
    vtkTemplateMacro(vtkExponentialDecayCompute
      (this, 
       this->FirstIteration, this->DecayFactor, this->LastUpdateTime, timevalue,
       output, arrays, Nc, static_cast<VTK_TT *>(0)));
  default:
    vtkErrorMacro(<< "Execute: Unknown ScalarType");
    return 0;
  }

  return output;
}
//----------------------------------------------------------------------------
vtkDataSet *vtkNeuronSpikeFilter::DecayDataSet(vtkDataSet *in1, double timevalue)
{
  vtkDataSet *output = in1->NewInstance();
  output->CopyStructure(in1);
  output->GetPointData()->ShallowCopy(in1->GetPointData());
  output->GetCellData()->ShallowCopy(in1->GetCellData());

  //
  // Compute spatial difference if the dataset is a vtkPointSet
  //
  std::vector<vtkDataArray*> arrays;
  vtkDataArray *outarray;
  //
  //
  // Loop over all pointdata 
  //
  for (int s=0; s<in1->GetPointData()->GetNumberOfArrays(); ++s) {
    arrays.clear();
    //
    // On some data, the scalar arrays are consistent but ordered
    // differently on each time step, so we will fetch them by name if
    // possible.
    //
    vtkDataArray *dataarray = in1->GetPointData()->GetArray(s);
    char *scalarname = dataarray->GetName();
    if (this->GetPointArrayStatus(scalarname)) {
      arrays.push_back(dataarray);
      vtkDataArray *dataarray = LastPointData->GetArray((std::string(this->ArrayNamePrefix)+std::string(scalarname)).c_str());
      arrays.push_back(dataarray); // NULL is OK
      outarray = this->DecayDataArray(timevalue, &arrays[0], arrays[0]->GetNumberOfTuples());
      this->FirstIteration = false;
      output->GetPointData()->AddArray(outarray);
      outarray->FastDelete();
    }
  }


//  if (in1->GetInformation()->Has(vtkDataObject::DATA_GEOMETRY_UNMODIFIED()))
//  {
//    output->GetInformation()->Set(vtkDataObject::DATA_GEOMETRY_UNMODIFIED(),1);
//  }

  return output;
}
*/
/*
//----------------------------------------------------------------------------
int vtkNeuronSpikeFilter::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
const char* vtkNeuronSpikeFilter::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkNeuronSpikeFilter::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkNeuronSpikeFilter::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name))
  {
    //    this->MeshParamsModifiedTime.Modified();
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
*/