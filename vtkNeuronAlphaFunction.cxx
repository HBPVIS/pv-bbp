/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkNeuronAlphaFunction.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkNeuronAlphaFunction.h"

#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkSmartPointer.h"
#include <cmath>
#include <string>
//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkNeuronAlphaFunction);
//----------------------------------------------------------------------------
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
//----------------------------------------------------------------------------
#define FloatOrDouble(F, D, index) F ? F[index] : D[index]
#define FloatOrDoubleorDefault(F, D, def, index) F ? F[index] : (D ? D[index] : def)
#define FloatOrDoubleSet(F, D) ((F!=NULL) || (D!=NULL))
//----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
vtkNeuronAlphaFunction::vtkNeuronAlphaFunction()
{
  this->DifferentialBlendFactor  = 1.0;
  this->HyperPolarizedVoltage    = -85.0;
  this->DePolarizedVoltage       = -50.0;
  this->RestingPotentialVoltage  = -65.0;
  this->PeakDifferentialVoltage  =  10.0;
  this->VoltageTransparencyMode  = 0;
  //
  this->Array1Name = NULL;
  this->Array2Name = NULL;
  this->Array3Name = NULL;
  this->Array4Name = NULL;
}
//-----------------------------------------------------------------------------
vtkNeuronAlphaFunction::~vtkNeuronAlphaFunction()
{
  delete []this->Array1Name;
  delete []this->Array2Name;
  delete []this->Array3Name;
  delete []this->Array4Name;
}
//-----------------------------------------------------------------------------
int vtkNeuronAlphaFunction::RequestData(
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

  //
  double *a1D, *a2D, *a3D, *a4D;
  float  *a1F, *a2F, *a3F, *a4F;
  // Get name of data array containing alpha
  vtkDataArray *a1 = this->Array1Name ? input->GetPointData()->GetArray(this->Array1Name) : NULL;
  FloatOrDoubleArrayPointer(a1, a1F, a1D);

  vtkDataArray *a2 = this->Array2Name ? input->GetPointData()->GetArray(this->Array2Name) : NULL;
  FloatOrDoubleArrayPointer(a2, a2F, a2D);

  vtkDataArray *a3 = this->Array3Name ? input->GetPointData()->GetArray(this->Array3Name) : NULL;
  FloatOrDoubleArrayPointer(a3, a3F, a3D);

  vtkDataArray *a4 = this->Array4Name ? input->GetPointData()->GetArray(this->Array4Name) : NULL;
  FloatOrDoubleArrayPointer(a4, a4F, a4D);

  vtkSmartPointer<vtkFloatArray> newAlpha = vtkSmartPointer<vtkFloatArray>::New();
  newAlpha->SetName("NeuronAlpha");

  vtkIdType N = a1->GetNumberOfTuples();
  float *newdata = newAlpha->WritePointer(0, N);

  // HyperPolarized transparent
  if (this->VoltageTransparencyMode==0) {
    for (vtkIdType i=0; i<N; i++) {
      double value = FloatOrDouble(a1F,a1D,i);
      double alpha = (value-this->HyperPolarizedVoltage)/(this->DePolarizedVoltage-this->HyperPolarizedVoltage);
      newdata[i] = (alpha<0 ? 0 : (alpha>1.0 ? 1.0 : alpha));
    }
  }

  // RestingPotential transparent
  if (this->VoltageTransparencyMode==1) {
    for (vtkIdType i=0; i<N; i++) {
      double value = FloatOrDouble(a1F,a1D,i);
      if (value>this->RestingPotentialVoltage) {
        double alpha = (value-this->RestingPotentialVoltage)/(this->DePolarizedVoltage-this->RestingPotentialVoltage);
        newdata[i] = (alpha<0 ? 0 : (alpha>1.0 ? 1.0 : alpha));
      }
      else {
        double alpha = (value-this->RestingPotentialVoltage)/(this->HyperPolarizedVoltage-this->RestingPotentialVoltage);
        newdata[i] = (alpha<0 ? 0 : (alpha>1.0 ? 1.0 : alpha));
      }
    }
  }

  if (a2 && a3) {
    for (vtkIdType i=0; i<N; i++) {
      double   rtopacity = FloatOrDouble(a2F,a2D,i);
      double dvdtopacity = FloatOrDouble(a3F,a3D,i);
      //
      double alpha = this->DifferentialBlendFactor*dvdtopacity + (1.0 - this->DifferentialBlendFactor)*rtopacity;
      newdata[i] = (alpha<0 ? 0 : (alpha>1.0 ? 1.0 : alpha));
    }
  }

//  for (vtkIdType i=0; i<N; i++) {
//    double value = a1->GetTuple1(i);
//    newdata[i] = 1.0 - ((1.0-value)*this->DifferentialBlendFactor);
//  }

  output->GetPointData()->AddArray(newAlpha);
  return 1;
}
