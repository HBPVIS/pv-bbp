// Visual studio debug settings for paths specifi to this module
//
// PATH=D:\build\paraview\bin\Debug;C:\Program Files\hdf5-vfd-1.8.9\bin;%PATH%
// PV_PLUGIN_PATH=D:\build\buildyard\ParaBBP\bin\Debug
//

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
#include "vtkTransform.h"

//
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDataArraySelection.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
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
#include "vtkOutlineSource.h"
#include "vtkAppendPolyData.h"
#include "vtkPolyDataNormals.h"
//
#include <vtksys/SystemTools.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vector>
#include <deque>
//
#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCharArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"
#include "vtkExtentTranslator.h"
#include "vtkBoundingBox.h"
#include "vtkIdListCollection.h"
//
#include "vtkDummyController.h"
//
#include <functional>
#include <algorithm>
#include <numeric>
#include <iterator>

// BBP-SDK Morphology Reader
#include "BBP/IO/File/Parsers/Morphology_HDF5_File_Parser.h"

// Header of the Reader
#include "vtkCircuitReader.h"

// BBP-SDK

#include <BBP/common.h>
#include "BBP/Microcircuit/Morphology.h"
#include "BBP/Microcircuit/Experiment.h"
#include "BBP/Microcircuit/Targets/Targets.h"
#include "BBP/Microcircuit/Targets/Cell_Target.h"
#include <BBP/Microcircuit/Containers/Neurons.h>

// Voxelization
#include "BBP/Voxelization/voxelization.h"
#include "BBP/VtkDebugging/visualization.h"


using namespace bbp;
//----------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkCircuitReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkCircuitReader);
//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkCircuitReader, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkCircuitReader::vtkCircuitReader()
{
  this->FileName = NULL;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(2);
  //
  this->NumberOfTimeSteps               = 0;
  this->TimeStep                        = 0;
  this->ActualTimeStep                  = 0;
  this->TimeStepTolerance               = 1E-6;
  this->FileName                        = NULL;
  this->DefaultTarget                   = NULL;
  this->UpdatePiece                     = 0;
  this->UpdateNumPieces                 = 0;
  this->IntegerTimeStepValues           = 0;
  this->MeshType                        = 0;
  this->GenerateNormalVectors           = 0;
  this->Random                          = 1;
  this->MaximumNumberOfNeurons          = 25;
  //
  this->PointDataArraySelection         = vtkDataArraySelection::New();
  this->TargetsSelection                = vtkDataArraySelection::New();
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
  this->SIL                      = vtkSmartPointer<vtkMutableDirectedGraph>::New();
  this->SILUpdateStamp           = 0;
}
//----------------------------------------------------------------------------
vtkCircuitReader::~vtkCircuitReader()
{
  this->SIL      = NULL;
  delete []this->FileName;
  this->PointDataArraySelection->Delete();
  this->TargetsSelection->Delete();
  this->SetController(NULL);
}
//----------------------------------------------------------------------------
int vtkCircuitReader::FillOutputPortInformation( int port, vtkInformation* info )
{
  if (port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
    return 1;
  }

  if (port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
    return 1;
  }

  return 0;
} 
//----------------------------------------------------------------------------
int vtkCircuitReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{

  // if there is no blue config supplied yet, exit quietly
  if (!this->FileName) {
    return 1;
  }

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  bool NeedToReadInformation = (FileModifiedTime>FileOpenedTime);

  if (NeedToReadInformation) {

    // ----------------------------------------------------------------------
    // Set parameters
    // ----------------------------------------------------------------------

    std::string blueconfig = this->FileName;

    // -------------------------------------------------------------------   
    // Create BBP-SDK Experiment and Microcircuit to access to the neurons.
    // -------------------------------------------------------------------   
    experiment.open(blueconfig);
    bbp::Microcircuit& microcircuit = experiment.microcircuit();

    this->NumberOfTimeSteps = 1;
    this->TimeStepValues.assign(this->NumberOfTimeSteps, 0.0);
    int validTimes = 0;
    for (int i=0; i<this->NumberOfTimeSteps; ++i)
    {
    }

    if (this->NumberOfTimeSteps==0) {
      vtkErrorMacro(<<"No time steps in data");
      return 0;
    }

    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
      &this->TimeStepValues[0],
      static_cast<int>(this->TimeStepValues.size()));
    double timeRange[2];
    timeRange[0] = this->TimeStepValues.front();
    timeRange[1] = this->TimeStepValues.back();
    if (this->TimeStepValues.size()>1)
    {
      this->TimeStepTolerance = 0.01*(this->TimeStepValues[1]-this->TimeStepValues[0]);
    }
    else
    {
      this->TimeStepTolerance = 1E-3;
    }
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

    this->BuildSIL();
    outInfo->Set(vtkDataObject::SIL(), this->GetSIL());

  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkCircuitReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{

  // open the HDF5 file
  if (!this->FileName) {
    vtkErrorMacro(<< "A BlueConfig FileName must be specified.");
    return 0;
  }

  // get the info object
  vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

  // get the ouptut
  vtkPolyData *output0 = vtkPolyData::SafeDownCast(outInfo0->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output1 = vtkPolyData::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));

  // default target?
  std::string target_name = this->DefaultTarget ? this->DefaultTarget : this->TargetsSelection->GetArrayName(0);
  //
  int N = this->TargetsSelection->GetNumberOfArrays();
  for (int i=0; i<N; i++) {
    const char *name = this->TargetsSelection->GetArrayName(i);
    if (this->TargetsSelection->ArrayIsEnabled(name)) {
      target_name = name;
      break;
      //std::cout << "Target selected : " << name << std::endl;
    }
  }

  bbp::Microcircuit& microcircuit = experiment.microcircuit();
  bbp::Target target;
  try {
    target = experiment.user_targets().get_target(target_name);
  }
  catch (...) {
    target = experiment.targets().get_target(target_name);
  }
  microcircuit.load(target, bbp::NEURONS | bbp::MORPHOLOGIES | bbp::MESHES);
  bbp::Neurons & neurons = microcircuit.neurons(); 
  //
  std::cout << "Neuron count for target : " << target_name.c_str() << " is " << neurons.size() << std::endl;

  //
  // Allocate VTK arrays
  //
  vtkSmartPointer<vtkPoints>       points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIntArray>   neuronId = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkFloatArray> nvectors = vtkSmartPointer<vtkFloatArray>::New();
  neuronId->SetName("NeuronId");
  nvectors->SetName("Normals");
  nvectors->SetNumberOfComponents(3);
  vtkSmartPointer<vtkPolyData>       normpoly = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>       normpoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> normtriangles = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
  normals->SetSplitting(0);
  normals->SetComputeCellNormals(0);
  normals->SetNonManifoldTraversal(0);
  //
  vtkSmartPointer<vtkTransform>     transform = vtkSmartPointer<vtkTransform>::New();
  vtkSmartPointer<vtkMatrix4x4>        matrix = vtkSmartPointer<vtkMatrix4x4>::New();

  vtkSmartPointer<vtkPoints>       pointsM = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray>     linesM = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPointData> pointdataM = output1->GetPointData();

  // -------------------------------------------------------------------   
  // iterate over neurons
  // -------------------------------------------------------------------   

#define USE_BBP_TRANSFORM
#undef  USE_VTK_TRANSFORM

  vtkIdType Ncount = 0;
  vtkIdType maxP = 0;
  vtkIdType maxF = 0;
  
  // count up the vertices and faces before allocating memory
  for (bbp::Neurons::iterator ni = neurons.begin(); ni != neurons.end(); ++ni,++Ncount) {
    if (this->MaximumNumberOfNeurons>0 && Ncount>=this->MaximumNumberOfNeurons) break;
    //
    const bbp::Mesh *sdk_mesh = &ni->morphology().mesh();
    maxP += sdk_mesh->vertex_count();
    maxF += sdk_mesh->triangle_count();
  }
  points->GetData()->Resize(maxP);
  points->SetNumberOfPoints(maxP);
  //    
  neuronId->Resize(maxP);
  neuronId->SetNumberOfTuples(maxP);
  //
  if (this->GenerateNormalVectors) {
    nvectors->Resize(maxP);
    nvectors->SetNumberOfTuples(maxP);
  }
  //
  vtkIdType *cells = triangles->WritePointer(maxF, 4*(maxF));

  //
  //
  //
  // Track the current insertion position for vertices
  vtkIdType insertN = 0;
  // each new neuron counts vertices from 0, we must increment by a growing offset with each neuron
  vtkIdType offsetN = 0;
  // Track insertion location for Cell vertex index Ids
  vtkIdType insertC = 0;
  //

  Ncount = 0;
  for (bbp::Neurons::iterator ni = neurons.begin(); ni != neurons.end(); ++ni,++Ncount) {
    
    // only load the maximum requested to save memory
    if (this->MaximumNumberOfNeurons>0 && Ncount>=this->MaximumNumberOfNeurons) break;
    //
    const bbp::Morphology* sdk_morph    = &ni->morphology();
    const bbp::Mesh*       sdk_mesh     = &ni->morphology().mesh();
    const bbp::Sections*   sdk_sections = &ni->neurites();

    const Transform_3D<bbp::Micron> &bbp_transform = ni->global_transform();
    this->CreateDatasetFromMorphology(sdk_morph, pointsM, linesM, pointdataM, bbp_transform);


    //
    bbp::Vertex_Index           vertexCount = sdk_mesh->vertex_count();
    bbp::Triangle_Index           faceCount = sdk_mesh->triangle_count();
    bbp::Triangle_Index         stripLength = sdk_mesh->triangle_strip_length();
    const Vector_3D<bbp::Micron> *vertices2 = sdk_mesh->vertices().pointer(); 
    //
    std::cout << "Neuron : " << std::distance(neurons.begin(),ni) << std::endl;
    std::cout << "vertex count : " << vertexCount << std::endl;
    std::cout << "face count : " << faceCount << std::endl;

    //
    // we must create a polydata object for the neuron in order to run the normal filter on it
    // better to do one neuron at a time than all of them in one go at the end.
    //
    vtkIdType *ncells = NULL;
    if (this->GenerateNormalVectors) {
      normpoints->GetData()->Resize(vertexCount);
      normpoints->SetNumberOfPoints(vertexCount);
      // resizing cell array doesn't work when new size is smaller - so create a new one
      normtriangles = vtkSmartPointer<vtkCellArray>::New();
      ncells = normtriangles->WritePointer(faceCount, 4*(faceCount));
      normpoly->SetPoints(normpoints);
      normpoly->SetPolys(normtriangles);
    }

#ifdef USE_VTK_TRANSFORM
    Degree_Angle y_rotation = ni->orientation().rotation;
    transform->PostMultiply();
    transform->RotateY(y_rotation);
    transform->Translate(ni->position().x(), ni->position().y(), ni->position().z());
    transform->Update();
    for (bbp::Vertex_Index v=0 ; v<vertexCount; ++v) {
      float newPoint[3];
      transform->TransformPoint(vertices2[v].vector(),newPoint); 
      points->SetPoint(insertN, newPoint);
      if (this->GenerateNormalVectors) {
        normpoints->setPoint(v, newPoint);
      }
#else
    for (bbp::Vertex_Index v=0 ; v<vertexCount; ++v) {
      bbp::Vector_3D<bbp::Micron> newPoint = bbp_transform*vertices2[v];
      points->SetPoint(insertN, newPoint.vector());
      if (this->GenerateNormalVectors) {
        normpoints->SetPoint(v, newPoint.vector());
      }
#endif
      neuronId->SetValue(insertN,Ncount);
      insertN++;
    }

    //
    // Generate triangles in cell array for each face
    // 
    const bbp::Vertex_Index* faces2 = sdk_mesh->triangles().pointer();
    for (bbp::Triangle_Index t=0; t<faceCount; ++t)
    {
      cells[insertC*4 + 0] = 3;
      cells[insertC*4 + 1] = offsetN + faces2[t*3 + 0];
      cells[insertC*4 + 2] = offsetN + faces2[t*3 + 1];
      cells[insertC*4 + 3] = offsetN + faces2[t*3 + 2];
      if (this->GenerateNormalVectors) {
        // no offset here as we only do one neuron at a time
        ncells[t*4 + 0] = 3;
        ncells[t*4 + 1] = faces2[t*3 + 0];
        ncells[t*4 + 2] = faces2[t*3 + 1];
        ncells[t*4 + 3] = faces2[t*3 + 2];
      }
      insertC++;
    }

    if (this->GenerateNormalVectors) {
      normals->Modified();
      normals->SetInput(normpoly);
      normals->Update();
      insertN = offsetN;
      vtkFloatArray *nvecs = vtkFloatArray::SafeDownCast(normals->GetOutput()->GetPointData()->GetArray("Normals"));
      if (nvecs->GetNumberOfTuples()!=vertexCount) {
          std::cout << "These numbers don't add up " << std::endl;
      }
      for (bbp::Vertex_Index v=0 ; v<vertexCount; ++v) {
        nvectors->SetTuple(insertN++, nvecs->GetTuple(v));
      }
    }
    offsetN = insertN;
  }
  //
  output0->SetPoints(points);
  output0->SetPolys(triangles);
  output0->GetPointData()->AddArray(neuronId);

  output1->SetPoints(pointsM);
  output1->SetLines(linesM);

  if (this->GenerateNormalVectors) {
    output0->GetPointData()->SetNormals(nvectors);
  }

  return 1;
}
//-----------------------------------------------------------------------------
void vtkCircuitReader::CreateDatasetFromMorphology(const bbp::Morphology *morph, vtkPoints *points, vtkCellArray *lines, vtkFieldData *field, const Transform_3D<bbp::Micron> &transform)
{
  // get the ouptut
  vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();

  // Load Morphology informations using the BBP-SDK parser
  Morphology_Repair_Stage repair_stage = REPAIRED_MORPHOLOGY;

  Morphology_Dataset dataset = std::move(morph->operator Morphology_Dataset());

  const Morphology_Point_ID  *section_start_points = dataset.section_start_points();
  const Section_Type                *section_types = dataset.section_types();
  const Section_ID                *section_parents = dataset.section_parent_sections();
  const Section_ID                   section_count = dataset.number_of_sections();
  const Vector_3D<bbp::Micron>    *point_positions = dataset.point_positions();
  const bbp::Micron               *point_diameters = dataset.point_diameters();
  const Morphology_Point_ID            point_count = dataset.point_count();

  std::cout << "section_count = " << section_count << std::endl;
  std::cout << "point_count = " << point_count << std::endl;

  vtkIdType offsetN = points->GetNumberOfPoints();
  vtkIdType maxP = offsetN + point_count;
  points->GetData()->Resize(maxP);
  points->SetNumberOfPoints(maxP);
  for (bbp::Morphology_Point_ID i=0; i<point_count; i++) {
    
    bbp::Vector_3D<bbp::Micron> newPoint = transform*point_positions[i];
    points->SetPoint(offsetN + i, newPoint.vector());
  }

  // Create the lines (cells) for the morphology
  vtkIdType offsetC = lines->GetNumberOfCells();
  vtkIdType maxL = offsetC;
  for (bbp::Section_ID sec_id=0; sec_id<section_count; sec_id++) {
    bbp::Morphology_Point_ID last_sec_pt = (sec_id < section_count - 1) ? section_start_points[sec_id+1] - 1: point_count - 1;
    for (bbp::Morphology_Point_ID pt_id=section_start_points[sec_id]; pt_id<last_sec_pt; pt_id++) { maxL++; }
  }
  vtkIdType *cells = lines->WritePointer(maxL, 3*(maxL));

  vtkIdType insertL = offsetC;
  for (bbp::Section_ID sec_id=0; sec_id<section_count; sec_id++) {
    bbp::Morphology_Point_ID last_sec_pt = (sec_id < section_count - 1) ? section_start_points[sec_id+1] - 1: point_count - 1;
    for (bbp::Morphology_Point_ID pt_id=section_start_points[sec_id]; pt_id<last_sec_pt; pt_id++)
    {
      cells[insertL*3 + 0] = 2;
      cells[insertL*3 + 1] = pt_id + offsetN;
      cells[insertL*3 + 2] = pt_id + offsetN + 1;
      insertL++;
    }
  }

  // Varying dendrite radius using the information of the morphology points
  vtkSmartPointer<vtkDoubleArray> dendriteRadius = vtkDoubleArray::SafeDownCast(field->GetArray("DendriteRadius"));
  if (!dendriteRadius) {
    dendriteRadius = vtkSmartPointer<vtkDoubleArray>::New();
    dendriteRadius->SetName("DendriteRadius");
    field->AddArray(dendriteRadius);
  }
  vtkSmartPointer<vtkUnsignedCharArray> sectionTypes = vtkUnsignedCharArray::SafeDownCast(field->GetArray("SectionTypes"));
  if (!sectionTypes) {
    sectionTypes = vtkSmartPointer<vtkUnsignedCharArray>::New();
    sectionTypes->SetName("SectionTypes");
    field->AddArray(sectionTypes);
  }
  vtkSmartPointer<vtkIntArray> sectionIds = vtkIntArray::SafeDownCast(field->GetArray("SectionIds"));
  if (!sectionIds) {
    sectionIds = vtkSmartPointer<vtkIntArray>::New();
    sectionIds->SetName("SectionIds");
    field->AddArray(sectionIds);
  }
  dendriteRadius->Resize(maxP);
  dendriteRadius->SetNumberOfTuples(maxP);
  sectionTypes->Resize(maxP);
  sectionTypes->SetNumberOfTuples(maxP);
  sectionIds->Resize(maxP);
  sectionIds->SetNumberOfTuples(maxP);
  //
  for (bbp::Morphology_Point_ID i=0; i<point_count ; i++)
    dendriteRadius->SetValue(offsetN + i, point_diameters[i]/2);
  // Overwrite the radius for the soma point (radius=0)
  for (bbp::Morphology_Point_ID i=0; i<section_start_points[1]-1; i++)
    dendriteRadius->SetValue(offsetN + i, 0.5);

  // Set a flag for each section type and Id
  bbp::Section_ID sec_id = 0;
  for (bbp::Morphology_Point_ID i=0; i<point_count; i++) {
    if (sec_id < section_count -1)
      if (i == section_start_points[sec_id+1])
        sec_id++;

    sectionIds->SetValue(offsetN + i,sec_id);
    bbp::Section_Type type = section_types[sec_id];
    switch(type) { // don't need the switch here, but useful to see flags
      case bbp::SOMA:
        sectionTypes->SetValue(offsetN + i,type);
        break;
      case bbp::AXON:
        sectionTypes->SetValue(offsetN + i,type);
        break;
      case bbp::DENDRITE:
        sectionTypes->SetValue(offsetN + i,type);
        break;
      case bbp::APICAL_DENDRITE:
        sectionTypes->SetValue(offsetN + i,type);
        break;
      default:
        sectionTypes->SetValue(offsetN + i,type);
    }
  }
}
//-----------------------------------------------------------------------------
void vtkCircuitReader::BuildSIL()
{
  // Initialize the SIL, dump all previous information.
  this->SIL->Initialize();

  vtkSmartPointer<vtkVariantArray> childEdge = vtkSmartPointer<vtkVariantArray>::New();
  childEdge->InsertNextValue(0);

  vtkSmartPointer<vtkVariantArray> crossEdge = vtkSmartPointer<vtkVariantArray>::New();
  crossEdge->InsertNextValue(0);

  // CrossEdge is an edge linking hierarchies.
  vtkSmartPointer<vtkUnsignedCharArray> crossEdgesArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  crossEdgesArray->SetName("CrossEdges");
  this->SIL->GetEdgeData()->AddArray(crossEdgesArray);

  std::deque<std::string> names;
  std::deque<vtkIdType> crossEdgeIds; // ids for edges between trees.
  int cc;

  // Now build the hierarchy.
  vtkIdType rootId = this->SIL->AddVertex();
  names.push_back("SIL");
  
  // Add the Heirarchy subtree.
  vtkIdType targetsRoot = this->SIL->AddChild(rootId, childEdge);
  names.push_back("Targets");

  // Get default targets for the microcircuit.
  const bbp::Targets &default_targets = this->experiment.targets();
  for (bbp::Targets::const_iterator ti=default_targets.begin(); ti!=default_targets.end(); ++ti) {
    std::string name = (*ti).name();
    vtkIdType childBlock = this->SIL->AddChild(targetsRoot, childEdge);
    names.push_back(name.c_str());
    this->TargetsSelection->AddArray(name.c_str());
  }

  // Get user targets for the microcircuit.
  const bbp::Targets &user_targets = this->experiment.user_targets();
  for (bbp::Targets::const_iterator ti=user_targets.begin(); ti!=user_targets.end(); ++ti) {
    std::string name = (*ti).name();
    vtkIdType childBlock = this->SIL->AddChild(targetsRoot, childEdge);
    names.push_back(name.c_str());
    this->TargetsSelection->AddArray(name.c_str());
  }

  // This array is used to assign names to nodes.
  vtkSmartPointer<vtkStringArray> namesArray = vtkSmartPointer<vtkStringArray>::New();
  namesArray->SetName("Names");
  namesArray->SetNumberOfTuples(this->SIL->GetNumberOfVertices());
  this->SIL->GetVertexData()->AddArray(namesArray);

  std::deque<std::string>::iterator iter;
  for (cc=0, iter = names.begin(); iter != names.end(); ++iter, ++cc)
  {
    const char *name = (*iter).c_str();
    namesArray->SetValue(cc, name);
  }
}
//----------------------------------------------------------------------------
void vtkCircuitReader::SetFileName(char *filename)
{
  if (this->FileName == NULL && filename == NULL)
  {
    return;
  }
  if (this->FileName && filename && (!strcmp(this->FileName,filename)))
  {
    return;
  }
  delete [] this->FileName;
  this->FileName = NULL;

  if (filename)
  {
    this->FileName = vtksys::SystemTools::DuplicateString(filename);
    this->SetFileModified();
  }
  this->Modified();
}
//----------------------------------------------------------------------------
void vtkCircuitReader::SetFileModified()
{
  this->FileModifiedTime.Modified();
  this->Modified();
}
//-----------------------------------------------------------------------------
vtkSmartPointer<vtkMutableDirectedGraph> vtkCircuitReader::GetSIL()
{
  return this->SIL;
}
//----------------------------------------------------------------------------
int vtkCircuitReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
const char* vtkCircuitReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkCircuitReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkCircuitReader::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name))
  {
    if (status)
    {
      this->PointDataArraySelection->EnableArray(name);
    }
    else
    {
      this->PointDataArraySelection->DisableArray(name);
    }
    this->Modified();
  }
}
//----------------------------------------------------------------------------
int vtkCircuitReader::GetNumberOfTargets()
{
  return this->TargetsSelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
const char* vtkCircuitReader::GetTargetsName(int index)
{
  return this->TargetsSelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkCircuitReader::GetTargetsStatus(const char* name)
{
  return this->TargetsSelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkCircuitReader::SetTargetsStatus(const char* name, int status)
{
  if(status)
    {
    this->TargetsSelection->EnableArray(name);
    }
  else
    {
    this->TargetsSelection->DisableArray(name);
    }
}
//----------------------------------------------------------------------------
void vtkCircuitReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
    << (this->FileName ? this->FileName : "(none)") << "\n";  
}
//----------------------------------------------------------------------------
