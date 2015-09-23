/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCircuitReaderMesh.cxx,v $

  Copyright (c) Sebastien Lasserre, Blue Brain Project
  All rights reserved.

=========================================================================*/

#ifndef __vtkCircuitReaderMesh_h
#define __vtkCircuitReaderMesh_h 
 
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

//
// BBP-SDK
//
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
class VTK_EXPORT vtkCircuitReaderMesh : public vtkCircuitReaderBase
{
public:
  static vtkCircuitReaderMesh *New();
  vtkTypeMacro(vtkCircuitReaderMesh, vtkCircuitReaderBase);

  void SetMaximumNumberOfNeurons(int n) {
    if (this->MaximumNumberOfNeurons != n) { 
      this->MaximumNumberOfNeurons = n; 
      this->MeshParamsModifiedTime.Modified();
      this->Modified(); 
    } 
  }
  vtkGetMacro(MaximumNumberOfNeurons,int);

  void SetExportMorphologySkeleton(int n) {
    if (this->ExportMorphologySkeleton != n) { 
      this->ExportMorphologySkeleton = n; 
      this->MeshParamsModifiedTime.Modified();
      this->Modified(); 
    } 
  }
  vtkGetMacro(ExportMorphologySkeleton,int);
  vtkBooleanMacro(ExportMorphologySkeleton,int);

  void SetExportNeuronMesh(int n) {
    if (this->ExportNeuronMesh != n) { 
      this->ExportNeuronMesh = n; 
      this->MeshParamsModifiedTime.Modified();
      this->Modified(); 
    } 
  }
  vtkGetMacro(ExportNeuronMesh,int);
  vtkBooleanMacro(ExportNeuronMesh,int);

  // Description:
  // hyperpolarized (voltage near -85 mV)
  vtkSetMacro(HyperPolarizedVoltage,double);
  vtkGetMacro(HyperPolarizedVoltage,double);

  // Description:
  // depolarized (voltage near -50 mV)
  vtkSetMacro(DePolarizedVoltage,double);
  vtkGetMacro(DePolarizedVoltage,double);

  // Description:
  // resting potential (voltage near -65 mV)
  vtkSetMacro(RestingPotentialVoltage,double);
  vtkGetMacro(RestingPotentialVoltage,double);
  
protected:
   vtkCircuitReaderMesh();
  ~vtkCircuitReaderMesh();

  int  RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  // top level functions which generate meshes/scalars for the whole data
  void GenerateNeuronMesh(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  void GenerateMorphologySkeleton(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  void CreateReportScalars(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  // generate mesh/scalars on a per neuron basis
  void      AddOneNeuronToMesh(bbp::Neuron *neuron, const bbp::Mesh *mesh, vtkIdType Ncount, vtkPoints *points, vtkIdType *cells, vtkFieldData *field, vtkIdType &offsetN, vtkIdType &offsetC);
  void      AddOneNeuronToMorphologySkeleton(bbp::Neuron *neuron, vtkIdType Ncount, vtkPoints *points, vtkIdType *cells, vtkFieldData *field, vtkIdType &offsetN, vtkIdType &offsetC);
  vtkIdType AddReportScalarsToMorphologySkeleton(bbp::Neuron *neuron, vtkFloatArray *voltages, vtkIdType offsetN);
  vtkIdType AddReportScalarsToNeuronMesh(bbp::Neuron *neuron, vtkFloatArray *voltages, vtkIdType offsetN);

  int ExportMorphologySkeleton;
  int ExportNeuronMesh;
  int MaximumNumberOfNeurons;

  double HyperPolarizedVoltage;
  double DePolarizedVoltage;
  double RestingPotentialVoltage;

  vtkSmartPointer<vtkPolyData>               CachedNeuronMesh;
  vtkSmartPointer<vtkPolyData>               CachedMorphologySkeleton;

private:
  vtkCircuitReaderMesh(const vtkCircuitReaderMesh&); 
  void operator=(const vtkCircuitReaderMesh&); 
};
 
#endif

