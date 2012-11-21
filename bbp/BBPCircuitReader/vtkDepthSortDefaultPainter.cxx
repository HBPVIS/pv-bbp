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
#include "vtkScalarsToColorsPainter.h"
#include "vtkGarbageCollector.h"
#include "vtkClipPlanesPainter.h"
#include "vtkPointsPainter.h"
#include "vtkStandardPolyDataPainter.h"
#include "vtkDepthSortPainter.h"
#include "vtkMapper.h"

vtkStandardNewMacro(vtkDepthSortDefaultPainter)

vtkCxxSetObjectMacro(vtkDepthSortDefaultPainter, DepthSortPainter, vtkDepthSortPainter)

vtkDepthSortDefaultPainter::vtkDepthSortDefaultPainter()
{
  this->DepthSortPainter = vtkDepthSortPainter::New();
}

vtkDepthSortDefaultPainter::~vtkDepthSortDefaultPainter()
{
  this->SetDepthSortPainter(NULL);
}

void vtkDepthSortDefaultPainter::BuildPainterChain()
{
  // build the superclass painter chain
  this->Superclass::BuildPainterChain();

  // insert the point sprite painter chain :
  // ScalarsToColorsPainter -> DepthSortPainter -> ...
  this->DepthSortPainter->SetDelegatePainter( this->ScalarsToColorsPainter->GetDelegatePainter());
  this->ScalarsToColorsPainter->SetDelegatePainter( this->DepthSortPainter);
}

void vtkDepthSortDefaultPainter::ReportReferences(vtkGarbageCollector *collector)
{
  this->Superclass::ReportReferences(collector);

  vtkGarbageCollectorReport(collector, this->DepthSortPainter,
      "DepthSortPainter");
}

void vtkDepthSortDefaultPainter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "DepthSortPainter: "
      << this->DepthSortPainter << endl;
}
