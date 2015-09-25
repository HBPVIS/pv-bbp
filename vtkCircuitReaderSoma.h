/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCircuitReaderSoma.cxx,v $

  Copyright (c) Sebastien Lasserre, Blue Brain Project
  All rights reserved.

=========================================================================*/

#ifndef __vtkCircuitReaderSoma_h
#define __vtkCircuitReaderSoma_h 
 
#include "vtkCircuitReaderBase.h"
#include "vtkSmartPointer.h"
#include "vtkUnsignedIntArray.h"
#define VTK_WRAPPING_CXX
#include "vtkClientServerStream.h"
#include <string>
#include <vector>


class vtkDataArraySelection;
class vtkMultiProcessController;
class vtkFloatArray;
class vtkMutableDirectedGraph;
class vtkMeshPartitionFilter;
class vtkUnstructuredGrid;
class vtkBoundsExtentTranslator;
// BBP-SDK

//#include "BBP/BBP.h"
#include "BBP/Targets/Target.h"
#include "BBP/Targets/Targets.h"
#include "BBP/Datasets/compartmentReportFrame.h"
#include "BBP/Readers/compartmentReportReader.h"
#include "BBP/Experiment.h"

//
#undef vtkGetMacro
#undef vtkSetMacro
//
#define vtkGetMacro(name,type) \
virtual type Get##name () { \
  return this->name; \
  }

#define vtkSetMacro(name,type) \
virtual void Set##name (type _arg) \
  { \
  if (this->name != _arg) \
    { \
    this->name = _arg; \
    this->Modified(); \
    } \
  }
//
class VTK_EXPORT vtkCircuitReaderSoma : public vtkCircuitReaderBase
{
public:
  static vtkCircuitReaderSoma *New();
  //vtkTypeMacro(vtkCircuitReaderSoma,vtkPolyDataAlgorithm);
  vtkTypeMacro(vtkCircuitReaderSoma, vtkCircuitReaderBase);

protected:
   vtkCircuitReaderSoma();
  ~vtkCircuitReaderSoma();

  // store infor about somas
  typedef std::tuple<double, vmml::Vector3f, int, unsigned char> soma_info;
  typedef std::map<uint32_t, soma_info> soma_map_type;

  int  RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  // top level functions which generate meshes/scalars for the whole data
  void GenerateSomaPoints(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  vtkSmartPointer<vtkPolyData>  CachedNeuronSoma;
  soma_map_type                 soma_map;

private:
  vtkCircuitReaderSoma(const vtkCircuitReaderSoma&); 
  void operator=(const vtkCircuitReaderSoma&); 
};
 
#endif

