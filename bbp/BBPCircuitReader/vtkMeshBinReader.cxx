
#include "vtkCellArray.h" 
#include "vtkTriangle.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"

// BBP-SDK Mesh Reader
#include "BBP/IO/File/Parsers/Mesh_Binary_File_Parser.h"

// Header of the Reader
#include "vtkMeshBinReader.h"

using namespace bbp;

vtkStandardNewMacro(vtkMeshBinReader);
 
vtkMeshBinReader::vtkMeshBinReader()
{
  this->FileName = NULL;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}
 
vtkMeshBinReader::~vtkMeshBinReader()
{
 
}
 
int vtkMeshBinReader::FillOutputPortInformation( int port, vtkInformation* info )
{
  if ( port == 0 )
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
 
    return 1;
  }
 
  return 0;
}
 
int vtkMeshBinReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
 
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
 
  // get the ouptut
   vtkPolyData *output = vtkPolyData::SafeDownCast(
            outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // open the binary mesh file
  if (!this->FileName)
  {
    vtkErrorMacro(<< "A FileName must be specified.");
    return 0;
  }

  vtkDebugMacro(<<"Reading file");

   Filename bin_file = this->FileName;
   bool vertexFlag = true;
   bool faceFlag = true;
   bool stripFlag = false;
   bool mappingFlag = false;
   Vertex_Index vertexCount;
   Triangle_Index faceCount;
   Triangle_Index stripLength;
   Vector_3D_Micron_Array vertices;
   Section_ID_Array vertex_sections;
   Float_Array vertex_relative_distances;
   Vertex_Index_Array faces;
   Vertex_Index_Array strips;

   Mesh_Binary_File_Parser binary_parser;
   binary_parser.read_mesh
   (
      bin_file,
      vertexFlag,
      faceFlag,
      stripFlag,
      mappingFlag,
      vertexCount,
      faceCount,
      stripLength,
      vertices,
      vertex_sections,
      vertex_relative_distances,
      faces,
      strips
   );
  
  // Set VTK points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  
  const Vector_3D<Micron>* vertices2 = vertices.get(); 
  std::cout << "vertex count : " << vertexCount << std::endl;
  for ( Vertex_Index v = 0 ; v < vertexCount ; ++v )
  {
     points->InsertNextPoint(vertices2[v].x(),
                             vertices2[v].y(),
                             vertices2[v].z());
  }
  
  // Set VTK triangles
  vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

  const Vertex_Index* faces2 = faces.get();
  std::cout  << faceCount << std::endl;
  for ( Triangle_Index t = 0 ; t < faceCount ; ++t )
  {
     vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
     triangle->GetPointIds()->SetId ( 0, faces[t*3]);
     triangle->GetPointIds()->SetId ( 1, faces[t*3+1]);
     triangle->GetPointIds()->SetId ( 2, faces[t*3+2]);
     triangles->InsertNextCell ( triangle );
  }
  
  // Set vtkPolyData
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetPolys(triangles);

  //output = polydata; //doesn't work
  output->ShallowCopy(polydata);
 
  return 1;
}
 
 
//----------------------------------------------------------------------------
void vtkMeshBinReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
 
  os << indent << "File Name: " 
      << (this->FileName ? this->FileName : "(none)") << "\n";  
}
