/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkH5MCellReader.h
  Revision of last commit : $Rev: 884 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2010-04-06 12:03:55 +0200 #$

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
#include "vtkH5MCellReader.h"
//
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDataArraySelection.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
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
#include "vtkCellArray.h"
#include "vtkOutlineSource.h"
#include "vtkAppendPolyData.h"
#include "vtkBoundingBox.h"
//
#include <vtksys/SystemTools.hxx>
#include <vtksys/RegularExpression.hxx>
//
#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCharArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"
#include "vtkExtentTranslator.h"
#include "vtkIdListCollection.h"
//
#include "vtkDummyController.h"
//
#include <algorithm>
#include <functional>
#include <numeric>
//
//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkH5MCellReader, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
/*!
  \ingroup h5part_utility

  This function can be used to query the Type of a dataset
  It is not used by the core H5Part library but is useful when
  reading generic data from the file.
  An example of usage would be (H5Tequal(datatype,H5T_NATIVE_FLOAT))
  any NATIVE type can be used to test.

  \return  \c an hdf5 handle to the native type of the data
*/
static hid_t H5PartGetNativeDatasetType(hid_t loc, const char *name)
{
  hid_t dataset, datatype, datatypen;
#if (!H5_USE_16_API && ((H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))))
  dataset=H5Dopen(loc, name, H5P_DEFAULT);
#else
  dataset=H5Dopen(loc, name);
#endif
  datatype  = H5Dget_type(dataset);
  datatypen = H5Tget_native_type(datatype, H5T_DIR_DEFAULT);
  H5Tclose(datatype);
  H5Dclose(dataset);
  return datatypen;
}

//----------------------------------------------------------------------------
hid_t H5PartGetDiskShape(hid_t file, hid_t dataset)
{
  hid_t space = H5Dget_space(dataset);
/*
  if (H5PartHasView(f))
    {
    int r;
    hsize_t stride, count;
    hsize_t range[2];
    // so, is this selection inclusive or exclusive? 
    range[0]=f->viewstart;
    range[1]=f->viewend;
    count = range[1]-range[0]; // to be inclusive
    stride=1;
    // now we select a subset 
    if (f->diskshape>0)
      {
      r=H5Sselect_hyperslab(f->diskshape,H5S_SELECT_SET,
        range, // only using first element
        &stride,&count,NULL);
      }
    // now we select a subset
    r=H5Sselect_hyperslab(space,H5S_SELECT_SET,
      range,&stride,&count,NULL);
    if (r<0)
      {
      fprintf(stderr,"Abort: Selection Failed!\n");
      return space;
      }
    }
*/
  return space;
}
//----------------------------------------------------------------------------
//#define JB_DEBUG__
#ifdef JB_DEBUG__
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
vtkStandardNewMacro(vtkH5MCellReader);
//----------------------------------------------------------------------------
vtkH5MCellReader::vtkH5MCellReader()
{
  this->SetNumberOfInputPorts(0);
  //
  this->NumberOfTimeSteps               = 0;
  this->TimeStep                        = 0;
  this->ActualTimeStep                  = 0;
  this->TimeStepTolerance               = 1E-6;
  this->GenerateVertexCells             = 0;
  this->FileName                        = NULL;
  this->H5FileId                        = NULL;
  this->UpdatePiece                     = 0;
  this->UpdateNumPieces                 = 0;
  this->TimeOutOfRange                  = 0;
  this->MaskOutOfTimeRangeOutput        = 0;
  this->IntegerTimeStepValues           = 0;
  this->PointDataArraySelection         = vtkDataArraySelection::New();
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
}
//----------------------------------------------------------------------------
vtkH5MCellReader::~vtkH5MCellReader()
{
  this->CloseFile();
  delete [] this->FileName;
  this->FileName = NULL;
 
  this->PointDataArraySelection->Delete();
  this->PointDataArraySelection = 0;

  this->SetController(NULL);
}
//----------------------------------------------------------------------------
void vtkH5MCellReader::SetFileName(char *filename)
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
void vtkH5MCellReader::SetFileModified()
{
  this->FileModifiedTime.Modified();
  this->Modified();
}
//----------------------------------------------------------------------------
void vtkH5MCellReader::CloseFile()
{
  
  if (this->H5FileId != NULL)
    {
//    H5PartCloseFile(this->H5FileId);
    this->H5FileId = NULL;
    }
}
//----------------------------------------------------------------------------
void vtkH5MCellReader::CloseFileIntermediate()
{
  H5Fclose(this->H5FileId);
  this->H5FileId = 0;
}
//----------------------------------------------------------------------------
int vtkH5MCellReader::OpenFile()
{
  if (!this->FileName)
    {
    vtkErrorMacro(<<"FileName must be specified.");
    return 0;
    }

  if (FileModifiedTime>FileOpenedTime)
    {
    this->CloseFile();
    }

  if (!this->H5FileId)
    {
    this->H5FileId = H5Fopen(this->FileName, H5F_ACC_RDONLY,H5P_DEFAULT);
    this->FileOpenedTime.Modified();
    }

  if (!this->H5FileId)
    {
    vtkErrorMacro(<< "Initialize: Could not open file " << this->FileName);
    return 0;
    }

  return 1;
}
//----------------------------------------------------------------------------
int GetVTKDataType(int datatype)
{
  if (H5Tequal(datatype,H5T_NATIVE_FLOAT))
    {
    return VTK_FLOAT;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_DOUBLE))
    {
    return VTK_DOUBLE;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_SCHAR))
    {
    return VTK_CHAR;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_UCHAR))
    {
    return VTK_UNSIGNED_CHAR;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_SHORT))
    {
    return VTK_SHORT;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_USHORT))
    {
    return VTK_UNSIGNED_SHORT;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_INT))
    {
    return VTK_INT;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_UINT))
    {
    return VTK_UNSIGNED_INT;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_LONG))
    {
    return VTK_LONG;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_ULONG))
    {
    return VTK_UNSIGNED_LONG;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_LLONG))
    {
    return VTK_LONG_LONG;
    }
  else if (H5Tequal(datatype,H5T_NATIVE_ULLONG))
    {
    return VTK_UNSIGNED_LONG_LONG;
    }
  return VTK_VOID;
}

//----------------------------------------------------------------------------
vtkSmartPointer<vtkDataArray> vtkH5MCellReader::ReadDataSetIntoDataArray(char *name, vtkIdType s, vtkIdType e)
{
  hid_t datatype = H5PartGetNativeDatasetType(H5FileId,name);
  int vtk_datatype = GetVTKDataType(datatype);
  if (vtk_datatype == VTK_VOID) {
    H5Tclose(datatype);
    vtkErrorMacro("An unexpected data type was encountered");
    return 0;
  }
  //
  // Get the array size
  //
  herr_t herr, r;
  hid_t dataset_id = H5Dopen(this->H5FileId, name, H5P_DEFAULT);
  hid_t  diskshape = H5Dget_space ( dataset_id );
  hsize_t dims[2], maxdims[2];
  herr_t err = H5Sget_simple_extent_dims(diskshape, dims, maxdims);
  int numDims = H5Sget_simple_extent_ndims(diskshape); 
  vtkIdType Nt = dims[0];
  vtkIdType Nc = numDims>1 ? dims[1] : 1;

  // if user has requested a piece
  if (s!=-1 && e!=-1) {
    Nt = (e-s);
    hsize_t  count1_mem[] = { Nt*Nc,  1 }; // the physical memory is really a flat array of size Nt*Nc 
    hsize_t  count2_mem[] = { Nt,    Nc }; // the memory space is Nt*Nc
    hsize_t  offset_mem[] = { 0,     0 };  // always read into meory starting at 0,0
    hsize_t  stride_mem[] = { 1,     1 };  // stride just 1,1
    hsize_t  count1_dsk[] = { Nt,   Nc };  // disk space is Nt,Nc
    hsize_t  offset_dsk[] = { s,     0 };  // offset from our start point
    hsize_t  stride_dsk[] = { 1,     1 };  // and stride 1,1
    r = H5Sselect_hyperslab(diskshape, H5S_SELECT_SET, offset_dsk, stride_dsk, count1_dsk, NULL);
  }
  vtkSmartPointer<vtkDataArray> dataarray;
  dataarray.TakeReference(vtkDataArray::CreateDataArray(vtk_datatype));
  dataarray->SetNumberOfComponents(Nc);
  dataarray->SetNumberOfTuples(Nt);
  dataarray->SetName(vtksys::SystemTools::GetFilenameName(name).c_str());
  herr = H5Dread(dataset_id, datatype, H5S_ALL, diskshape, H5P_DEFAULT, dataarray->GetVoidPointer(0));

  herr = H5Tclose(datatype);
  herr = H5Sclose ( diskshape );
  herr = H5Dclose ( dataset_id );

  return dataarray;
}
//----------------------------------------------------------------------------
int vtkH5MCellReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  bool NeedToReadInformation = (FileModifiedTime>FileOpenedTime || !this->H5FileId);

  if (NeedToReadInformation)
    {
    if (!this->OpenFile())
      {
      return 0;
      }
    this->offsets = this->ReadDataSetIntoDataArray("/molecules/glutamate/0/positions/offsets", -1, -1);

    this->NumberOfTimeSteps = this->offsets->GetNumberOfTuples();
/*
    for (int i=0; i<nds; i++)
      {
      // return 0 for no, 1,2,3,4,5 etc for index (1 based offset)
      H5PartGetDatasetName(this->H5FileId, i, name, 512);
      this->PointDataArraySelection->AddArray(name);
      }
*/
    this->TimeStepValues.assign(this->NumberOfTimeSteps, 0.0);

    if (this->NumberOfTimeSteps==0)
    {
      vtkErrorMacro(<<"No time steps in data");
      return 0;
    }

    // if TIME information was either not present ot not consistent, then
    // set something so that consumers of this data can iterate sensibly
    if (1 || this->IntegerTimeStepValues || (this->NumberOfTimeSteps>0))
      {
      for (int i=0; i<this->NumberOfTimeSteps; i++)
        {
        // insert read of Time array here
        this->TimeStepValues[i] = i;
        }
      }
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
      &this->TimeStepValues[0],
      static_cast<int>(this->TimeStepValues.size()));
    double timeRange[2];
    timeRange[0] = this->TimeStepValues.front();
    timeRange[1] = this->TimeStepValues.back();
    if (this->TimeStepValues.size()>1)
      {
      this->TimeStepTolerance = 0.01*(this->TimeStepValues[1]-this->TimeStepValues[0]);
      }
    else
      {
      this->TimeStepTolerance = 1E-3;
      }
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

  }

  this->CloseFileIntermediate();

  return 1;
}
//----------------------------------------------------------------------------
template <class T1, class T2>
void CopyIntoTuple(int offset, vtkDataArray *source, vtkDataArray *dest)
{
  vtkIdType N = source->GetNumberOfTuples();
  T1 *sptr = static_cast<T1*>(source->GetVoidPointer(0));
  T2 *dptr = static_cast<T2*>(dest->WriteVoidPointer(0,N)) + offset;
  for (vtkIdType i=0; i<N; ++i) {
    *dptr = *sptr++;
    dptr += 3;
  }
}
//----------------------------------------------------------------------------
/*
std::pair<double, double> GetClosest(std::vector<double> &sortedlist, const double& val) const
{
  std::vector<double>::const_iterator it = std::lower_bound(sortedlist.begin(), sortedlist.end(), val);
  if (it == sortedlist.end())        return std::make_pair(sortedlist.back(), sortedlist.back());
  else if (it == sortedlist.begin()) return std::make_pair(sortedlist.front(), sortedlist.front());
  else return std::make_pair(*(it - 1), *(it));
}
*/
class H5MCellToleranceCheck: public std::binary_function<double, double, bool>
{
public:
  H5MCellToleranceCheck(double tol) { this->tolerance = tol; }
  double tolerance;
  //
    result_type operator()(first_argument_type a, second_argument_type b) const
    {
      bool result = (fabs(a-b)<=(this->tolerance));
      return (result_type)result;
    }
};
//----------------------------------------------------------------------------
#if (!H5_USE_16_API && ((H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))))
  #define h_params ,H5P_DEFAULT
#else
  #define h_params
#endif
//----------------------------------------------------------------------------
int vtkH5MCellReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData     *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  //
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  //
  if (this->TimeStepValues.size()==0) return 0;
  //
//  int N = this->PointDataArraySelection->GetNumberOfArrays();
  //
  // Get the TimeStep Requested from the information if present
  //
  this->TimeOutOfRange = 0;
  this->ActualTimeStep = this->TimeStep;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
    {
    double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    this->ActualTimeStep = std::find_if(
      this->TimeStepValues.begin(), this->TimeStepValues.end(),
      std::bind2nd( H5MCellToleranceCheck( 
          this->IntegerTimeStepValues ? 0.5 : this->TimeStepTolerance ), requestedTimeValue ))
      - this->TimeStepValues.begin();
    //
    if (requestedTimeValue<this->TimeStepValues.front() || requestedTimeValue>this->TimeStepValues.back())
      {
      this->TimeOutOfRange = 1;
      }
    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);
    }
  else
    {
    double timevalue[1];
    unsigned int index = this->ActualTimeStep;
    if (index<this->TimeStepValues.size())
      {
      timevalue[0] = this->TimeStepValues[index];
      }
    else
      {
      timevalue[0] = this->TimeStepValues[0];
      }
    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), timevalue[0]);
    }

  if (this->TimeOutOfRange && this->MaskOutOfTimeRangeOutput)
    {
    // don't do anything, just return success
    return 1;
    }

  // open the file if not already done
  if (!this->OpenFile())
    {
    return 0;
    }


  // Setup arrays for reading data
  int startval = this->offsets->GetTuple1(this->TimeStep);
  int   endval = this->offsets->GetTuple1(this->TimeStep+1);
  vtkSmartPointer<vtkPoints>    points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkDataArray> coords = this->ReadDataSetIntoDataArray("/molecules/glutamate/0/positions/data", startval, endval);

  //
  // generate cells
  //
  vtkIdType Nt = coords->GetNumberOfTuples();
  points->SetNumberOfPoints(Nt);
  if (1 || this->GenerateVertexCells)
    {
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType *cells = vertices->WritePointer(Nt, 2*Nt);
    for (vtkIdType i=0; i<Nt; ++i)
      {
      cells[2*i] = 1;
      cells[2*i+1] = i;
      }
    output->SetVerts(vertices);
    }
  //
  //
  //
  points->SetData(coords);
  output->SetPoints(points);
  //
  //
  // only subclasses actually close the file.
  //
  this->CloseFileIntermediate();

  return 1;
}
//----------------------------------------------------------------------------
unsigned int mylog2(unsigned int val) {
  unsigned int ret = -1;
  while (val != 0) {
    val >>= 1;
    ret++;
  }
  return ret;
}
//----------------------------------------------------------------------------
const char* vtkH5MCellReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkH5MCellReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkH5MCellReader::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name))
    {
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
void vtkH5MCellReader::Enable(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}
//----------------------------------------------------------------------------
void vtkH5MCellReader::Disable(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}
//----------------------------------------------------------------------------
void vtkH5MCellReader::EnableAll()
{
  this->PointDataArraySelection->EnableAllArrays();
}
//----------------------------------------------------------------------------
void vtkH5MCellReader::DisableAll()
{
  this->PointDataArraySelection->DisableAllArrays();
}
//----------------------------------------------------------------------------
int vtkH5MCellReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
void vtkH5MCellReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " <<
    (this->FileName ? this->FileName : "(none)") << "\n";

  os << indent << "NumberOfSteps: " <<  this->NumberOfTimeSteps << "\n";
}
