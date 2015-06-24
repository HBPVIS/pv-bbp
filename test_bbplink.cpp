// Visual studio debug settings for paths specific to this module
// //
// // PATH=D:\build\paraview-3.98\bin\Debug;C:\Program Files\hdf5-vfd-1.8.9\bin;%PATH%
// // PV_PLUGIN_PATH=D:\build\buildyard\ParaBBP\bin\Debug
// // _NO_DEBUG_HEAP=1
// // working directory : D:\build\paraview-3.98\bin\Debug
// //
//
//



#include "vtkObjectFactory.h"
#include "vtkCellArray.h" 
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
//
#include "vtkDataArraySelection.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataArray.h"
//
#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkLongLongArray.h"
#include "vtkUnsignedLongLongArray.h"
#include "vtkIntArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkVariantArray.h"
#include "vtkStringArray.h"
#include "vtkCellArray.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkExtentTranslator.h"
//
#include "vtkTransform.h"
#include "vtkPolyDataNormals.h"
//
#include "vtkDummyController.h"
//
#include "vtkPKdTree.h"
#include "vtkBoundsExtentTranslator.h"
#include "vtkMeshPartitionFilter.h"
//
#include <vtksys/SystemTools.hxx>
//
#include <vector>
#include <deque>
#include <functional>
#include <algorithm>
#include <numeric>
#include <iterator>
//
// BBP-SDK Morphology Reader
#include "BBP/BBP.h"
//#include "BBP/Microcircuit/Readers/Microcircuit_Reader.h"
//#include "BBP/Microcircuit/Readers/Mesh_Reader.h"
//#include "BBP/Microcircuit/Readers/compartmentReportReader.h"
//#include "BBP/Microcircuit/Mappings/Compartment_Report_Mapping.h"
//#include "BBP/Microcircuit/Soma.h"
//#include "BBP/Microcircuit/Mesh.h"
//#include "BBP/Microcircuit/Datasets/Morphology_Dataset.h"

// Header of this Reader
#include "vtkCircuitReader.h"
//



int main(int argc, char **argv) {
  vtkSmartPointer<vtkCircuitReader> reader = vtkSmartPointer<vtkCircuitReader>::New();
}


