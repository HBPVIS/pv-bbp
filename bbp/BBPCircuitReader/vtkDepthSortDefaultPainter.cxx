/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDepthSortDefaultPainter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkDepthSortDefaultPainter
// .SECTION Thanks
// <verbatim>
//
//  This file is part of the DepthSorts plugin developed and contributed by
//
//  Copyright (c) CSCS - Swiss National Supercomputing Centre
//                EDF - Electricite de France
//
//  John Biddiscombe, Ugo Varetto (CSCS)
//  Stephane Ploix (EDF)
//
// </verbatim>

#include "vtkDepthSortDefaultPainter.h"

#include "vtkObjectFactory.h"
#include "vtkGarbageCollector.h"
#include "vtkMapper.h"
#include "vtkClipPlanesPainter.h"
#include "vtkPointsPainter.h"
#include "vtkStandardPolyDataPainter.h"
#include "vtkDepthSortPainter.h"
#include "vtkScalarsToColorsPainter.h"
#include "vtkTwoScalarsToColorsPainter.h"
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkDepthSortDefaultPainter)
//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkDepthSortDefaultPainter, DepthSortPainter, vtkDepthSortPainter)
vtkCxxSetObjectMacro(vtkDepthSortDefaultPainter, TwoScalarsToColorsPainter, vtkTwoScalarsToColorsPainter)
//----------------------------------------------------------------------------
vtkDepthSortDefaultPainter::vtkDepthSortDefaultPainter()
{
  this->DepthSortPainter          = vtkDepthSortPainter::New();
  this->TwoScalarsToColorsPainter = vtkTwoScalarsToColorsPainter::New();
}
//----------------------------------------------------------------------------
vtkDepthSortDefaultPainter::~vtkDepthSortDefaultPainter()
{
  this->SetDepthSortPainter(NULL);
  this->SetTwoScalarsToColorsPainter(NULL);
}
//----------------------------------------------------------------------------
void vtkDepthSortDefaultPainter::BuildPainterChain()
{
  // build the superclass painter chain
  this->Superclass::BuildPainterChain();

  // insert our painters into the current painter chain :
  // ... -> ScalarsToColorsPainter -> TwoScalarsToColorsPainter -> DepthSortPainter -> ...
  this->DepthSortPainter->SetDelegatePainter(this->ScalarsToColorsPainter->GetDelegatePainter());
  this->TwoScalarsToColorsPainter->SetDelegatePainter(this->DepthSortPainter);
  this->ScalarsToColorsPainter->SetDelegatePainter(this->TwoScalarsToColorsPainter);
}
//----------------------------------------------------------------------------
void vtkDepthSortDefaultPainter::ReportReferences(vtkGarbageCollector *collector)
{
  this->Superclass::ReportReferences(collector);

  vtkGarbageCollectorReport(collector, this->DepthSortPainter,
      "DepthSortPainter");
}
//----------------------------------------------------------------------------
void vtkDepthSortDefaultPainter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "DepthSortPainter: "
      << this->DepthSortPainter << endl;
}
//----------------------------------------------------------------------------
