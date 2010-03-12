#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkOutlineFilter.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPointData.h>
#include <vtkPiecewiseFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkVolumeTextureMapper3D.h>
#include <vtkVolume.h>
//#include <.h>
//#include <.h>
//#include <.h>
//#include <.h>

#include <bbp.h>

using namespace bbp;

int main (int argc, char *argv[])
{

  //parse command line arguments
  if(argc != 2)
  {
  vtkstd::cout << "Required arguments: Filename" << vtkstd::endl;
  exit(-1);
  }
    
  vtkstd::string filename = argv[1];

  Timer timer;
  Pipeline pipeline;

    const URI blue_config_file = filename;
    Micron  center_x(boost::lexical_cast<Micron>("275.0")),
            center_y(boost::lexical_cast<Micron>("600.0")),
            center_z(boost::lexical_cast<Micron>("275.0")),
            width(boost::lexical_cast<Micron>("55.0")),
            height(boost::lexical_cast<Micron>("120.0")),
            depth(boost::lexical_cast<Micron>("55.0")),
            resolution(boost::lexical_cast<Micron>("1.0"));
    // load experiment and microcircuit data
    Experiment experiment;
    experiment.open(blue_config_file);
    Microcircuit & microcircuit = experiment.microcircuit();

    Cell_Target cell_target = experiment.cell_target("MiniColumn_0");
#ifndef NDEBUG
    cell_target.print();
#endif
    timer.start("Loading microcircuit");
    microcircuit.load(cell_target, NEURONS | MORPHOLOGIES);
    timer.print();

    // define volume dimensions and resolution
    Box <Micron> bounding_box (center_x, center_y, center_z,  width, height, depth);

    // define the section_type
    const char* section_type_string = std::string("all").c_str();
    int sec_type = 0;
    if (!strcmp(section_type_string,"all") || !strcmp(section_type_string,"soma"))
        sec_type |= SOMA;
    if (!strcmp(section_type_string,"all") || !strcmp(section_type_string,"axons"))
        sec_type |= AXON;
    if (!strcmp(section_type_string,"all") || !strcmp(section_type_string,"dendrites"))
        sec_type |= DENDRITE | APICAL_DENDRITE;

    // configure and run processing pipeline
    timer.start("Setting up pipeline");
    bool clean_cut = true;
    Segment_Voxel_Mapper volume_mapping(resolution, bounding_box, clean_cut, Section_Type(sec_type));
    Segment_Voxel_Density density_volume;
    //Volume_NRRD_File_Stream_Writer nrrd_file_writer(std::string("output"));
    volume_mapping.input_ptr()->reset(experiment.microcircuit_ptr());

    pipeline = volume_mapping + density_volume;// + nrrd_file_writer;
    timer.print();

    timer.start("Running pipeline");
    pipeline.mode(PROCESS_SINGLE_STEP_TRIGGER);
    std::cout << "Starting pipeline" << std::endl;
    pipeline.start();
    // Process the SDK-Pipeline 
    std::cout << "Processing..." << std::endl;
    pipeline.step();
    timer.print();
  
    // Build the VTK image
    timer.start("Generate vtkImageData");
    vtkSmartPointer < vtkImageData > imageData = vtkSmartPointer<vtkImageData>::New();
    //specify the size of the image data
    imageData->SetDimensions(width,height,depth);
    imageData->SetNumberOfScalarComponents(1);
    imageData->SetScalarTypeToUnsignedChar();

    int* dims = imageData->GetDimensions();

    cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << endl;

    cout << "Number of points: " << imageData->GetNumberOfPoints() << endl;
    cout << "Number of cells: " << imageData->GetNumberOfCells() << endl;

    //fill every entry of the image data with "2.0"
    //std::vector<unsigned char> elts = density_volume.output().get_elements();
    //unsigned int index = 0;
    //index = x + y * dims[1] + z * dims[0] * dims[1];
    for (int z = 0; z < dims[2]; z++)
      {
        for (int y = 0; y < dims[1]; y++)
        {
        for (int x = 0; x < dims[0]; x++)
          {
          unsigned char* pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(x,y,z));
          pixel[0] = density_volume.output()(x,y,z);
          }
        }
      }
    timer.print();


    timer.start("Deconstruction");
    std::cout << "Stopping pipeline." << std::endl;
    pipeline.stop();
    timer.print();

// print image data content
  for (int z = 0; z < dims[2]; z++)
    {
    for (int y = 0; y < dims[1]; y++)
      {
      for (int x = 0; x < dims[0]; x++)
        {
        unsigned char* pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(x,y,z));
        cout << pixel[0] << " ";
        }
      cout << endl;
      }
    cout << endl;
    }


  // bounding box  
  vtkSmartPointer<vtkOutlineFilter> outline = vtkSmartPointer<vtkOutlineFilter>::New();  
  outline->SetInput( imageData );  
  vtkSmartPointer<vtkPolyDataMapper> outlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();  
  outlineMapper->SetInput( outline->GetOutput() );  
  vtkSmartPointer<vtkActor> outlineActor = vtkSmartPointer<vtkActor>::New();
  outlineActor->SetMapper( outlineMapper ); 


// volume
vtkPiecewiseFunction* opacityTransferFunction = vtkPiecewiseFunction::New();
   opacityTransferFunction->AddPoint(0, 0.0);
   opacityTransferFunction->AddPoint(5, 1.0);


vtkVolumeProperty* volumeProperty = vtkVolumeProperty::New();
   volumeProperty->SetScalarOpacity(opacityTransferFunction);
   volumeProperty->IndependentComponentsOff();
   

vtkVolumeTextureMapper3D* volumeMapper = vtkVolumeTextureMapper3D::New();
   volumeMapper->SetInput(imageData);     

vtkVolume* volume = vtkVolume::New();
   volume->SetMapper(volumeMapper);
   volume->SetProperty(volumeProperty);

std::cout << "Volume done ..." << std::endl;

  // a renderer and render window
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
 
  // an interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
 
  // add the actors to the scene
  renderer->AddActor(outlineActor);
  renderer->AddVolume(volume);
  renderer->SetBackground(.2, .3, .4); 
 
  // render an image (lights and cameras are created automatically)
  renderWindow->Render();

  // change the interactor style
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = 
    vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(style); 

  // begin mouse interaction
  renderWindowInteractor->Start();
 
  return 0;
}
