#include <vtkSmartPointer.h>
#include <vtkImageData.h>
 
int main(int argc, char *argv[])
{
  //create an image data
  vtkSmartPointer<vtkImageData> imageData = 
      vtkSmartPointer<vtkImageData>::New();
 
  //specify the size of the image data
  imageData->SetDimensions(2,3,1);
  imageData->SetNumberOfScalarComponents(1);
  imageData->SetScalarTypeToDouble();
 
  int* dims = imageData->GetDimensions();
  // int dims[3]; // can't do this
 
  cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << endl;
 
  cout << "Number of points: " << imageData->GetNumberOfPoints() << endl;
  cout << "Number of cells: " << imageData->GetNumberOfCells() << endl;
 
  //fill every entry of the image data with "2.0"
  for (int z = 0; z < dims[2]; z++)
    {
    for (int y = 0; y < dims[1]; y++)
      {
      for (int x = 0; x < dims[0]; x++)
        {
        double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
        pixel[0] = 2.0;
        }
      }
    }
 
  //retrieve the entries from the image data and print them to the screen
  for (int z = 0; z < dims[2]; z++)
    {
    for (int y = 0; y < dims[1]; y++)
      {
      for (int x = 0; x < dims[0]; x++)
        {
        double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
        // do something with v
        cout << pixel[0] << " ";
        }
      cout << endl;
      }
    cout << endl;
    }
 
  return EXIT_SUCCESS;
}
