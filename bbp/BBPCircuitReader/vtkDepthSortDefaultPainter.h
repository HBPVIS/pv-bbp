/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDepthSortDefaultPainter.h

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
// .SECTION Description
// The vtkDepthSortDefaultPainter replaces the vtkScalarsToColorsPainter by a vtkTwoScalarsToColorsPainter
// and add a vtkDepthSortPainter in the painter chain.

#ifndef __vtkDepthSortDefaultPainter_h__
#define __vtkDepthSortDefaultPainter_h__

#include "vtkDefaultPainter.h"

class vtkDepthSortPainter;
class vtkTwoScalarsToColorsPainter;
class vtkDepthSortPolygonsPainter;

class VTK_EXPORT vtkDepthSortDefaultPainter : public vtkDefaultPainter
{
public :
  static vtkDepthSortDefaultPainter* New();
  vtkTypeMacro(vtkDepthSortDefaultPainter, vtkDefaultPainter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get/Set the DepthSortPainter. 
  // The depth sort painter reorders primive in distance from camera order
  void SetDepthSortPainter(vtkDepthSortPainter*);
  vtkGetObjectMacro(DepthSortPainter, vtkDepthSortPainter);

  // Description:
  // Get/Set the TwoScalarsToColorsPainter. 
  // The two scalars painter combines the usual RGB scalars with another array
  // for the alpha component
  void SetTwoScalarsToColorsPainter(vtkTwoScalarsToColorsPainter*);
  vtkGetObjectMacro(TwoScalarsToColorsPainter, vtkTwoScalarsToColorsPainter);

  // Description:
  // Get/Set the DepthSortPolygonsPainter. 
  // The DepthSortPolygonsPainter is a primitive painter which paints
  // using the order provided by the sorted depths generated in the
  // Depth Sort Painter which should be earlier in the painter chain
  void SetDepthSortPolygonsPainter(vtkDepthSortPolygonsPainter*);
  vtkGetObjectMacro(DepthSortPolygonsPainter, vtkDepthSortPolygonsPainter);

protected:
  // Description:
  // Setups the the painter chain.
  virtual void BuildPainterChain();

  // Description:
  // Take part in garbage collection.
  virtual void ReportReferences(vtkGarbageCollector *collector);

  vtkDepthSortPainter          *DepthSortPainter;
  vtkTwoScalarsToColorsPainter *TwoScalarsToColorsPainter;
  vtkDepthSortPolygonsPainter  *DepthSortPolygonsPainter;

  vtkDepthSortDefaultPainter();
  ~vtkDepthSortDefaultPainter();

private:
  vtkDepthSortDefaultPainter(const vtkDepthSortDefaultPainter&); // Not implemented.
  void operator=(const vtkDepthSortDefaultPainter&); // Not implemented.
};

#endif
