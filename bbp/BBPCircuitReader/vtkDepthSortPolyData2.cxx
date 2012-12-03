/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkDepthSortPolyData2.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkDepthSortPolyData2.h"

#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPolyData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"

#include <algorithm>
#include <functional>
#include <tuple>
typedef std::tuple<double, vtkIdType, vtkIdType> depthInfo;
typedef std::vector<depthInfo> depthList;

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkDepthSortPolyData2);
//-----------------------------------------------------------------------------
vtkDepthSortPolyData2::vtkDepthSortPolyData2()
{
  this->Direction       = VTK_DIRECTION_BACK_TO_FRONT;
  this->DepthSortMode   = VTK_SORT_FIRST_POINT;
  this->FastPolygonMode = 1;
  this->SortingList     = new depthList;
}
//-----------------------------------------------------------------------------
vtkDepthSortPolyData2::~vtkDepthSortPolyData2()
{
}
//-----------------------------------------------------------------------------
void FloatOrDoubleArrayPointer(vtkDataArray *dataarray, float *&F, double *&D) {
  if (dataarray && vtkFloatArray::SafeDownCast(dataarray)) {
    F = vtkFloatArray::SafeDownCast(dataarray)->GetPointer(0);
    D = NULL;
  }
  if (dataarray && vtkDoubleArray::SafeDownCast(dataarray)) {
    D = vtkDoubleArray::SafeDownCast(dataarray)->GetPointer(0);
    F = NULL;
  }
  //
  if (dataarray && !F && !D) {
    vtkGenericWarningMacro(<< dataarray->GetName() << "must be float or double");
  }
}
//-----------------------------------------------------------------------------
template<typename T> 
void CentreBoundsFromPtIds(vtkIdType *pts, vtkIdType npts, T *points, T result[3]) 
{
  T bounds[6];
  T *pt = &points[pts[0]*3];
  bounds[0] = bounds[1] = pt[0];
  bounds[2] = bounds[3] = pt[1];
  bounds[4] = bounds[5] = pt[2];
  for (vtkIdType i=1; i<npts; i++) {
    T *pt = &points[pts[0]*3];
    bounds[0] = std::min(pt[0],bounds[0]);
    bounds[1] = std::max(pt[0],bounds[1]);
    bounds[2] = std::min(pt[1],bounds[2]);
    bounds[3] = std::max(pt[1],bounds[3]);
    bounds[4] = std::min(pt[2],bounds[4]);
    bounds[5] = std::max(pt[2],bounds[5]);
  }
  result[0] = (bounds[0]+bounds[1])/2.0;
  result[1] = (bounds[2]+bounds[3])/2.0;
  result[2] = (bounds[4]+bounds[5])/2.0;
}
//-----------------------------------------------------------------------------
int vtkDepthSortPolyData2::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // if the user has not requested fast mode, default to the standard sort
  if (!this->FastPolygonMode) {
    return this->vtkDepthSortPolyData::RequestData(request, inputVector, outputVector);
  }

  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Compute the sort vector from camera position
  double vectorD[3], origin[3];
  if ( this->Camera == NULL) {
    vtkErrorMacro(<<"Need a camera to sort");
    return 0;
  }
  this->ComputeProjectionVector(vectorD, origin);
  float vectorF[3] = { vectorD[0], vectorD[1], vectorD[2] };
  //
  vtkCellArray *polys = input->GetPolys();
  vtkPoints   *points = input->GetPoints();
  vtkIdType  numCells = polys->GetNumberOfCells();
  //
  // Points are always float or double, optimize access by avoiding copies
  // we copy the pointer and use the data directly.
  //
  float     *pointsF = NULL;
  double    *pointsD = NULL;
  if (numCells>0) {
    FloatOrDoubleArrayPointer(points->GetData(), pointsF, pointsD);
  }

  // allocate space for depth values and sorted order
  depthList *ListToSort = static_cast<depthList*>(this->SortingList);
  ListToSort->resize(numCells);
  // traverse polygon list and compute depth
  polys->InitTraversal();
  vtkIdType *zerooffset = polys->GetPointer();
  for (vtkIdType cellId=0; cellId<numCells; cellId++) {
    vtkIdType *pts;
    vtkIdType  npts;
    // get pointer to point Ids for this cell
    polys->GetNextCell(npts, pts);

    if ( this->DepthSortMode == VTK_SORT_FIRST_POINT )
    {
      // set depth using float/double operation
      if (pointsF) {
        float *x = &pointsF[pts[0]*3];
        std::get<0>(ListToSort->operator[](cellId)) = vtkMath::Dot(x,vectorF);
      }
      else {
        double *x = &pointsD[pts[0]*3];
        std::get<0>(ListToSort->operator[](cellId)) = vtkMath::Dot(x,vectorD);
      }
    }
    else if ( this->DepthSortMode == VTK_SORT_BOUNDS_CENTER )
    {
      if (pointsF) {
        float x[3];
        CentreBoundsFromPtIds<float>(pts, npts, pointsF, x);
        std::get<0>(ListToSort->operator[](cellId)) = vtkMath::Dot(x,vectorF);
      }
      else {
        double x[3];
        CentreBoundsFromPtIds<double>(pts, npts, pointsD, x);
        std::get<0>(ListToSort->operator[](cellId)) = vtkMath::Dot(x,vectorD);
      }
    }
    else // VTK_SORT_PARAMETRIC_CENTER )
    {
    }

    // set the cell ID to this cell
    std::get<1>(ListToSort->operator[](cellId)) = cellId;
    // set the offset to the cell {N,ptIds}
    std::get<2>(ListToSort->operator[](cellId)) = static_cast<vtkIdType>(pts-zerooffset-1);
  }
  this->UpdateProgress(0.20);

  // Sort the depths
  std::sort(ListToSort->begin(), ListToSort->end(), std::greater<depthInfo>());
  this->UpdateProgress(0.60);

  //  outCD->CopyAllocate(inCD);
  //  output->Allocate(tmpInput,numCells);
  vtkSmartPointer<vtkIdTypeArray> DepthOrder = vtkSmartPointer<vtkIdTypeArray>::New();
  DepthOrder->SetName("DepthOrder");
  DepthOrder->SetNumberOfComponents(2);
  DepthOrder->SetNumberOfTuples(numCells);
  for (vtkIdType cellId=0; cellId<numCells; cellId++) {
    // set values 1,2,3 
    vtkIdType tupledata[2] = {
      // tuple consists of "sorted cell Id", "Offset into cellArray list for {N,PtIds}"
      std::get<1>(ListToSort->operator[](cellId)),
      std::get<2>(ListToSort->operator[](cellId)),
    };
    DepthOrder->SetTupleValue(cellId, tupledata);
  }
  this->UpdateProgress(0.90);

  // Points are left alone
  output->ShallowCopy(input);
  output->GetCellData()->AddArray(DepthOrder);

  return 1;
}
//-----------------------------------------------------------------------------
