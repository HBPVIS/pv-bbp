/*=========================================================================

  Project                 : Icarus
  Module                  : vtkNeuronSpikeTableSource.cxx

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

  This work has received funding from the European Community's Seventh
  Framework Programme (FP7/2007-2013) under grant agreement 225967 “NextMuSE”

=========================================================================*/
#include "vtkNeuronSpikeTableSource.h"
//
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTable.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
//
// For PARAVIEW_USE_MPI
#include "vtkPVConfig.h"
#ifdef PARAVIEW_USE_MPI
 #include "vtkMPI.h"
 #include "vtkMPIController.h"
 #include "vtkMPICommunicator.h"
#endif
// Otherwise
#include "vtkMultiProcessController.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkNeuronSpikeTableSource, NameStrings, vtkStringList);
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkNeuronSpikeTableSource);
//----------------------------------------------------------------------------
vtkNeuronSpikeTableSource::vtkNeuronSpikeTableSource()
{
  // we override the input and instead generate values directly
  this->SetNumberOfInputPorts(0);
  this->NameStrings          = NULL;
  this->NumberOfTimeSteps    = 0;
  this->LastMaxTime          = 0.0;
  this->LastMinTime          = 0.0;
  this->Values               = vtkSmartPointer<vtkPointData>::New();
  this->TimeData             = vtkSmartPointer<vtkDoubleArray>::New();
  this->TimeData->SetName("TimeData");
  this->SpikeData            = vtkSmartPointer<vtkDoubleArray>::New();
  this->SpikeData->SetName("Spike Frequency");
  this->BinResolution        = 0.1;
  this->MaximumTableSize     = 500;
  this->PreviousBinStartTime = 0.0;
//  this->NameStrings = vtkStringList::New();
//  this->NameStrings->AddString("Spike Frequency");
}

//----------------------------------------------------------------------------
vtkNeuronSpikeTableSource::~vtkNeuronSpikeTableSource()
{
//  this->NameStrings->Delete();
}

//----------------------------------------------------------------------------
void vtkNeuronSpikeTableSource::ClearSpikeData()
{
}

//----------------------------------------------------------------------------
void vtkNeuronSpikeTableSource::SetSpikeData(vtkIdType N, signed char Ids[])
{
//  vtkWarningMacro("SetSpikeData - Type 1 " << N);
}

//----------------------------------------------------------------------------
void vtkNeuronSpikeTableSource::SetSpikeData(vtkIdType N, vtkClientServerStreamDataArg<signed char> &temp0)
{
//  vtkWarningMacro("SetSpikeData - Type 2 " << N);
  //
  spike_info_type *begin = (spike_info_type *) (temp0.operator signed char *());
  this->spikelist.assign(begin, begin + N/sizeof(spike_info_type));
  //
  this->Modified();
  this->SpikeListModifiedTime.Modified();
}

//----------------------------------------------------------------------------
void vtkNeuronSpikeTableSource::SetDataModified()
{
  this->DataModifiedTime.Modified();
  this->Modified();
}

//----------------------------------------------------------------------------
int vtkNeuronSpikeTableSource::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  int result = vtkTableAlgorithm::RequestInformation(request, inputVector, outputVector);
  //
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //
  int numT = this->TimeData->GetNumberOfTuples();
  if (numT>0) {
    double *tvals = this->TimeData->GetPointer(0);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), tvals, numT);
  }

  return result;
}

//----------------------------------------------------------------------------
int vtkNeuronSpikeTableSource::ProcessSpikeData(vtkTable *outdata)
{
  double firsttime = VTK_DOUBLE_MAX;
  double lasttime = 0.0; // this->LastMaxTime;
  //
  for (const spike_info_type &it : spikelist) {
    double time = std::get<1>(it);
    firsttime=std::min(firsttime, time);
    lasttime=std::max(lasttime, time);
    //
  }
  if (firsttime>lasttime) {
    firsttime = 0.0;
    lasttime = 0.0;
  }
  if (lasttime>=this->LastMaxTime) {
    this->LastMaxTime = lasttime;
  }
   // time was reset to zero or something like that (loop demo)
  else {
    this->Flush();
    this->LastMaxTime = lasttime;
  }

  //
  this->LastMinTime = firsttime;

  double real_start = static_cast<int>((lasttime/this->BinResolution)+0.5) - this->MaximumTableSize;
  real_start = (real_start > 0) ? real_start : 0;
  // the startpoint of bin 0
  real_start = real_start*this->BinResolution;
  int num_bins = (lasttime-real_start) / this->BinResolution;
  //
  std::vector<int> bins(num_bins,0);
  for (const spike_info_type &it : spikelist) {
    double time = std::get<1>(it);
    int bin_num = (time-real_start) / this->BinResolution;
    if (bin_num>=0 && bin_num<num_bins) {
      bins[bin_num]++;
    }
  }
  // we have summed the current spikes, but we also want the ones we had previously
  for (int i=0; i<this->TimeData->GetNumberOfTuples(); ++i) {
    double time = this->TimeData->GetValue(i);
    int bin_num = (time-real_start) / this->BinResolution;
    if (bin_num>=0 && bin_num<num_bins) {
      int val = this->SpikeData->GetValue(i);
      bins[bin_num] += val;
    }
  }

  this->TimeData = vtkSmartPointer<vtkDoubleArray>::New();
  this->TimeData->SetName("TimeData");
  this->TimeData->SetNumberOfTuples(num_bins);
  for (int i=0; i<num_bins; ++i) {
    double tval = real_start + this->BinResolution*i;
    this->TimeData->SetValue(i, tval);
  }
  this->SpikeData = vtkSmartPointer<vtkDoubleArray>::New();
  this->SpikeData->SetName("Spike Frequency");
  this->SpikeData->SetNumberOfTuples(num_bins);
  for (int i=0; i<num_bins; ++i) {
    this->SpikeData->SetValue(i, bins[i]);
  }

  vtkIdType numT = this->TimeData->GetNumberOfTuples();

  // to stop annoying "can't plot with one point" error messages
  // we will always do the first point twice if there's only one
  if (numT==1 && this->TimeData->GetNumberOfTuples()>0) {
    this->TimeData->InsertNextTuple(this->TimeData->GetTuple(0));
    this->SpikeData->InsertNextTuple(this->SpikeData->GetTuple(0));
  }

  return 1;
}

//---------------------------------------------------------------------------
int vtkNeuronSpikeTableSource::RequestData(
  vtkInformation *vtkNotUsed(information),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkTable* output = vtkTable::GetData(outputVector);
  vtkInformation* outputInfo = outputVector->GetInformationObject(0);

  // did simulation write time into steering values?
  double currenttime = this->LastMaxTime;
/*
  if (this->DsmManager->GetSteeringValues("TimeValue", 1, &currenttime)==H5FD_DSM_SUCCESS) {
    // value set in currenttime
  }
  else 
  */
  if (outputInfo && outputInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
  {
    currenttime = outputInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  }

  this->ProcessSpikeData(output);

  output->AddColumn(this->TimeData);
  output->AddColumn(this->SpikeData);
  

/*
  //
  // if the user rewound the animation, the time will be wrong, reset
  //
  if (currenttime<this->LatestTime) {
    this->Flush();
  }
  this->LatestTime = currenttime;
*/

  return 1;
}
//---------------------------------------------------------------------------
void vtkNeuronSpikeTableSource::Flush()
{
  this->Values->Initialize();
  this->TimeData->Initialize();
  this->SpikeData->Initialize();
}

//---------------------------------------------------------------------------
void vtkNeuronSpikeTableSource::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
//-----------------------------------------------------------------------------
