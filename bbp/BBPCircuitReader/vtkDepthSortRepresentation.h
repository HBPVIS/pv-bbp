/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile$

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkDepthSortRepresentation
// .SECTION Description
// vtkDepthSortRepresentation is an extension for vtkGeometryRepresentation
// that renders point-sprites at all point locations.

#ifndef __vtkDepthSortRepresentation_h
#define __vtkDepthSortRepresentation_h

#include "vtkGeometryRepresentation.h"
#include "vtkSmartPointer.h"

class vtkDepthSortPainter;
class vtkDepthSortDefaultPainter;
class vtkMultiProcessController;
class vtkBoundsExtentTranslator;

class VTK_EXPORT vtkDepthSortRepresentation : public vtkGeometryRepresentation
{
public:
  static vtkDepthSortRepresentation* New();
  vtkTypeMacro(vtkDepthSortRepresentation, vtkGeometryRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // When on (default), the representation tells the view to use the
  // partitioning information from the input structured grid for ordered
  // compositing. When off we let the view build its own ordering and
  // redistribute data as needed.
  void SetUseDataParititions(bool);
  vtkGetMacro(UseDataParititions, bool);

  // Description:
  // By default this filter uses the global controller,
  // but this method can be used to set another instead.
  virtual void SetController(vtkMultiProcessController*);


protected:
  vtkDepthSortRepresentation();
  ~vtkDepthSortRepresentation();

  // Description:
  virtual int RequestData(vtkInformation*,
    vtkInformationVector**, vtkInformationVector*);

  virtual int ProcessViewRequest(
    vtkInformationRequestKey* request_type,
    vtkInformation* inInfo, vtkInformation* outInfo);


  vtkDepthSortDefaultPainter* DepthSortDefaultPainter;
  vtkDepthSortPainter* DepthSortPainter;
  vtkSmartPointer<vtkBoundsExtentTranslator> BoundsTranslator;

  int UseDataParititions;
  //
  vtkMultiProcessController *Controller;
  //
  double GlobalDataBounds[6];

private:
  vtkDepthSortRepresentation(const vtkDepthSortRepresentation&); // Not implemented
  void operator=(const vtkDepthSortRepresentation&); // Not implemented

};

#endif
