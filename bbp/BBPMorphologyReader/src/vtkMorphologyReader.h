/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkMorphologyReader.cxx,v $

  Copyright (c) Sebastien Lasserre, Blue Brain Project
  All rights reserved.

=========================================================================*/

#ifndef __vtkMorphologyReader_h
#define __vtkMorphologyReader_h 
 
#include "vtkPolyDataAlgorithm.h"
 
class VTK_IO_EXPORT vtkMorphologyReader : public vtkPolyDataAlgorithm 
{
public:
  static vtkMorphologyReader *New();
  //vtkTypeMacro(vtkMorphologyReader,vtkPolyDataAlgorithm);
  vtkTypeRevisionMacro(vtkMorphologyReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  // Description:
  // Specify file name of the morphology file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
 
protected:
  vtkMorphologyReader();
  ~vtkMorphologyReader();

  int FillOutputPortInformation( int port, vtkInformation* info );
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  char *FileName;
 
private:
  vtkMorphologyReader(const vtkMorphologyReader&); 
  void operator=(const vtkMorphologyReader&); 
};
 
#endif

