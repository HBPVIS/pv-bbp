/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkAlphaTransform.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef __vtkAlphaTransform_h
#define __vtkAlphaTransform_h

#include "vtkPointSetAlgorithm.h"


class VTK_EXPORT vtkAlphaTransform : public vtkPointSetAlgorithm
{
public:
  static vtkAlphaTransform *New();
  vtkTypeMacro(vtkAlphaTransform,vtkPointSetAlgorithm);

  // Description:
  // 
  vtkSetMacro(BlendFactor,double);
  vtkGetMacro(BlendFactor,double);


protected:
  vtkAlphaTransform();
  ~vtkAlphaTransform();

  double BlendFactor;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkAlphaTransform(const vtkAlphaTransform&);  // Not implemented.
  void operator=(const vtkAlphaTransform&);  // Not implemented.
};

#endif
