/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDepthSortPolygonsPainter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkDepthSortPolygonsPainter - this painter paints polygons.
// .SECTION Description
// This painter renders Polys in vtkPolyData. It can render the polys
// in any representation (VTK_POINTS, VTK_WIREFRAME, VTK_SURFACE).

#ifndef __vtkDepthSortPolygonsPainter_h
#define __vtkDepthSortPolygonsPainter_h

#include "vtkPolygonsPainter.h"

class VTK_EXPORT vtkDepthSortPolygonsPainter : public vtkPolygonsPainter
{
public:
  static vtkDepthSortPolygonsPainter* New();
  vtkTypeMacro(vtkDepthSortPolygonsPainter, vtkPolygonsPainter);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkDepthSortPolygonsPainter();
  ~vtkDepthSortPolygonsPainter();

  // Description:
  // The actual rendering happens here. This method is called only when
  // SupportedPrimitive is present in typeflags when Render() is invoked.
  virtual int RenderPrimitive(unsigned long flags, vtkDataArray* n,
    vtkUnsignedCharArray* c, vtkDataArray* t, vtkRenderer* ren);

private:
  vtkDepthSortPolygonsPainter(const vtkDepthSortPolygonsPainter&); // Not implemented.
  void operator=(const vtkDepthSortPolygonsPainter&); // Not implemented.

};


#endif

