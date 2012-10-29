/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: mesh_bin_reader.cxx,v $

  Copyright (c) Sebastien Lasserre, Blue Brain Project
  All rights reserved.

=========================================================================*/

#ifndef __vtkMeshBinReader_h
#define __vtkMeshBinReader_h 
 
#include "vtkPolyDataAlgorithm.h"
 
class VTK_EXPORT vtkMeshBinReader : public vtkPolyDataAlgorithm 
{
public:
  static vtkMeshBinReader *New();
  //vtkTypeMacro(vtkMeshBinReader,vtkPolyDataAlgorithm);
  vtkTypeRevisionMacro(vtkMeshBinReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  // Description:
  // Specify file name of the mesh file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
 
protected:
  vtkMeshBinReader();
  ~vtkMeshBinReader();

  int FillOutputPortInformation( int port, vtkInformation* info );
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  char *FileName;
 
private:
  vtkMeshBinReader(const vtkMeshBinReader&); // Not implemented
  void operator=(const vtkMeshBinReader&); // Not implemented
};
 
#endif

