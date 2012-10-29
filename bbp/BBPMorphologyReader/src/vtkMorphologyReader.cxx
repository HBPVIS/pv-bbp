
#include "vtkCellArray.h" 
#include "vtkTriangle.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyLine.h>
#include <vtkLine.h>
#include <vtkCellData.h>

// BBP-SDK Morphology Reader
#include "BBP/IO/File/Parsers/Morphology_HDF5_File_Parser.h"

// Header of the Reader
#include "vtkMorphologyReader.h"

using namespace bbp;

vtkCxxRevisionMacro(vtkMorphologyReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkMorphologyReader);
 
vtkMorphologyReader::vtkMorphologyReader()
{
  this->FileName = NULL;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}
 
vtkMorphologyReader::~vtkMorphologyReader()
{
 
}
 
int vtkMorphologyReader::FillOutputPortInformation( int port, vtkInformation* info )
{
  if ( port == 0 )
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
 
    return 1;
  }
 
  return 0;
}
 
int vtkMorphologyReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
 
    // get the info object
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
 
    // get the ouptut
    vtkPolyData *output = vtkPolyData::SafeDownCast(
            outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // open the HDF5 file
    if (!this->FileName)
    {
      vtkErrorMacro(<< "A FileName must be specified.");
      return 0;
    }

    vtkDebugMacro(<<"Reading file");

    // Load Morphology informations using the BBP-SDK parser
    Morphology_HDF5_File_Parser parser;
    Morphology_Repair_Stage repair_stage = REPAIRED_MORPHOLOGY;

    Morphology_Point_ID_Array section_start_points;
    Section_Type_Array        section_types;
    Section_ID_Array          section_parents;
    Section_ID                section_count;
    Vector_3D_Micron_Array    point_positions;
    Micron_Array              point_diameters;
    Morphology_Point_ID       point_count;

    parser.open(this->FileName, repair_stage);

    section_start_points    = parser.section_start_points();
    section_types           = parser.section_types();
    section_parents         = parser.parent_sections();
    section_count           = parser.number_of_sections();
    point_positions         = parser.point_positions();
    point_diameters         = parser.point_diameters();
    point_count             = parser.number_of_points();

    double vX, vY, vZ;

    std::cout << "section_count = " << section_count << std::endl;
    std::cout << "point_count = " << point_count << std::endl;
    

    // Create points for the morphology
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(bbp::Morphology_Point_ID i  = 0; i < point_count; i++)
    {
        vX = point_positions[i].x();
        vY = point_positions[i].y();
        vZ = point_positions[i].z();
        points->InsertPoint(i, vX, vY, vZ);
    }

    // Create the lines (cells) for the morphology
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
 
    for (bbp::Section_ID sec_id = 0; sec_id < section_count; sec_id++)
    {
        bbp::Morphology_Point_ID last_sec_pt = (sec_id < section_count - 1) ? 
            section_start_points[sec_id+1] - 1: point_count - 1;
        for (bbp::Morphology_Point_ID pt_id = section_start_points[sec_id];
                                      pt_id < last_sec_pt ; pt_id++)
        {
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            line->GetPointIds()->SetId(0,pt_id);
            line->GetPointIds()->SetId(1,pt_id+1);
            cells->InsertNextCell(line);
        }
    }

    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(cells);

    // Varying dendrite radius using the information of the morphology points
    vtkSmartPointer<vtkDoubleArray> dendriteRadius = vtkSmartPointer<vtkDoubleArray>::New();
    dendriteRadius->SetName("DendriteRadius");
    dendriteRadius->SetNumberOfTuples(point_count);
    for (bbp::Morphology_Point_ID i=0 ;i<point_count ; i++)
        dendriteRadius->SetTuple1(i,point_diameters[i]/2);
    // Overwrite the radius for the soma point (radius=0)
    for (bbp::Morphology_Point_ID i=0 ;i<section_start_points[1]-1 ; i++)
        dendriteRadius->SetTuple1(i,0.5);
    
    polyData->GetPointData()->AddArray(dendriteRadius);
    polyData->GetPointData()->SetActiveScalars("DendriteRadius");

    // RBG array for all the points.
    // I tried to attach a color to the lines but it is not working with the Tube filter
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetName("Colors");
    colors->SetNumberOfComponents(3);
    colors->SetNumberOfTuples(point_count);
    bbp::Section_ID sec_id = 0;
    for (bbp::Morphology_Point_ID i = 0; i < point_count ;i++)
    {
        if (sec_id < section_count -1)
            if (i == section_start_points[sec_id+1])
                sec_id++;
        bbp::Section_Type type = section_types[sec_id];
        switch(type)
        {
            case bbp::SOMA:
                colors->InsertTuple3(i,0,255,0);
                break;
            case bbp::AXON:
                colors->InsertTuple3(i,0,0,255);
                break;
            case bbp::DENDRITE:
                colors->InsertTuple3(i,255,100,50);
                break;
            case bbp::APICAL_DENDRITE:
                colors->InsertTuple3(i,255,0,0);
                break;
            default:
                colors->InsertTuple3(i,255,255,255);
        }
    }
    polyData->GetPointData()->AddArray(colors);
  
  
    output->ShallowCopy(polyData);

    return 1;
}
 
 
//----------------------------------------------------------------------------
void vtkMorphologyReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
 
  os << indent << "File Name: " 
      << (this->FileName ? this->FileName : "(none)") << "\n";  
}
