/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkAlphaTransform.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkAlphaTransform.h"

#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkSmartPointer.h"
#include <cmath>
#include <string>
//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkAlphaTransform);
//-----------------------------------------------------------------------------
vtkAlphaTransform::vtkAlphaTransform()
{
  this->BlendFactor = 1.0;
}
//-----------------------------------------------------------------------------
vtkAlphaTransform::~vtkAlphaTransform()
{
}
//-----------------------------------------------------------------------------
int vtkAlphaTransform::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and ouptut
  vtkPointSet *input = vtkPointSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPointSet *output = vtkPointSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // everything will go through except the array we modify
  output->ShallowCopy(input);

  // Get name of data array containing alpha
  vtkDataArray  *alpha = this->GetInputArrayToProcess(0, inputVector);
  if (!alpha || alpha->GetNumberOfTuples()<1)
    {
    return 1;
    }

  vtkSmartPointer<vtkFloatArray> newAlpha = vtkSmartPointer<vtkFloatArray>::New();
  std::string name = alpha->GetName() ? alpha->GetName()  : "";
  newAlpha->SetName(std::string(name+"_Alpha").c_str());

  vtkIdType N = alpha->GetNumberOfTuples();
  float *newdata = newAlpha->WritePointer(0, N);

  for (vtkIdType i=0; i<N; i++) {
    double value = alpha->GetTuple1(i);
    newdata[i] = 1.0 - ((1.0-value)*this->BlendFactor);
  }

  output->GetPointData()->AddArray(newAlpha);
  return 1;
}
