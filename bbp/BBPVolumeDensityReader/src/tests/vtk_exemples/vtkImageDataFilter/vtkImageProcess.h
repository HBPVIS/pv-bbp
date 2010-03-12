/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageProcess.cxx,v $

  Copyright (c) Sebastien Lasserre, Blue Brain Project
  All rights reserved.

=========================================================================*/

#ifndef __vtkImageProcess_h
#define __vtkImageProcess_h 
 
#include "vtkImageAlgorithm.h"

class VTK_IO_EXPORT vtkImageProcess : public vtkImageAlgorithm 
{
public:
  static vtkImageProcess *New();
  //vtkTypeMacro(vtkImageProcess,vtkImageAlgorithm);
  vtkTypeRevisionMacro(vtkImageProcess,vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
protected:
  vtkImageProcess();
  ~vtkImageProcess();

  int FillOutputPortInformation( int port, vtkInformation* info );
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkImageProcess(const vtkImageProcess&); 
  void operator=(const vtkImageProcess&);

};
 
#endif

