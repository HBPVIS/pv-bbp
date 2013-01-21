#ifdef _WIN32
#include <windows.h>
#endif

#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>
#include <thrust/sort.h>

#include <cuda_gl_interop.h>

#include "vtkScalarsToColors.h"
#include "vtkPistonDataObject.h"
#include "vtkPistonDataWrangling.h"
#include "vtkPistonScalarsColors.h"
#include "vtkPistonMinMax.h"
#include "vtkPistonReference.h"
#include "piston/piston_math.h"

#include "vtkgl.h"

#include <iostream>

using namespace std;

namespace vtkpiston {

bool AlmostEqualRelativeAndAbs(float A, float B,
            float maxDiff, float maxRelDiff)
  {
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    float diff = fabs(A - B);
    if (diff <= maxDiff)
        return true;

    A = fabs(A);
    B = fabs(B);
    float largest = (B > A) ? B : A;

    if (diff <= largest * maxRelDiff)
        return true;
    return false;
  }

template <typename ValueType>
struct color_map : thrust::unary_function<ValueType, float3>
{
    const ValueType min;
    const ValueType max;
    const int size;
    float *table;
    const int numberOfChanels;

    color_map(float *table, int arrSize, int noOfChanels,
      ValueType rMin, ValueType rMax) :
      min(rMin),
      max(rMax),
      size((arrSize / noOfChanels) - 1),
      table(table),
      numberOfChanels(noOfChanels)
      {
      }

    __host__ __device__
    float3 operator()(ValueType val)
    {
      int index = 0;
      if((max - min) > 0.0)
        {
        index = ( (val - min) / (max - min) ) * size;
        }

      if (index < 0) index = 0;
      if (index > size) index = size;
      index *= numberOfChanels;

      float3 color;
      if(numberOfChanels == 1)
        {
        color = make_float3(table[index], table[index], table[index]);
        }
      else if(numberOfChanels == 2)
        {
        color = make_float3(table[index], table[index + 1], 0.0f);
        }
      else if(numberOfChanels == 3)
        {
        color = make_float3(table[index], table[index + 1], table[index + 2]);
        }
      else
        {
        // Not supported
        }

      return color;
    }
};

//------------------------------------------------------------------------------
void CudaGLInit()
{
  cudaDeviceProp prop;
  int dev;

  // Fill it with zeros
  memset(&prop,0,sizeof(cudaDeviceProp));

  // Pick a GPU capable of 1.0 or better
  prop.major=1; prop.minor=0;
  cudaChooseDevice(&dev,&prop);

  // Set OpenGL device
  cudaError_t res = cudaGLSetGLDevice(dev);

  if (res != cudaSuccess)
    {
    cerr << "Set device failed  ... " << cudaGetErrorString(res) << endl;
    return;
    }
}

//------------------------------------------------------------------------------
void CudaRegisterBuffer(struct cudaGraphicsResource **vboResource,
                        GLuint vboBuffer)
{
  cudaError_t res =
    cudaGraphicsGLRegisterBuffer(vboResource, vboBuffer,
                                cudaGraphicsMapFlagsWriteDiscard);
  if (res != cudaSuccess)
  {
    cerr << "Register buffer failed ... " << cudaGetErrorString(res) << endl;
    return;
  }
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
struct distance_functor 
{
  float3 cameravector;

  // construct with a constant camera vector
  __host__ __device__ distance_functor(float3 &cam) : cameravector(cam) {}

  template <typename Tuple>
  __host__ __device__
  void operator()(Tuple t)
  {
    thrust::get<1>(t) = dot(thrust::get<0>(t), cameravector);  
  }
};
//------------------------------------------------------------------------------
struct celldistance_functor 
{
  const float *vertex_distances;
  
  // construct with a precomputed distance vector for every vertex
  __host__ __device__ celldistance_functor(float *v) : vertex_distances(v) {}
  
  template <typename Tuple>
  __host__ __device__
  void operator()(Tuple t)
  { 
    thrust::get<1>(t) = (vertex_distances[thrust::get<0>(t).x] + 
                         vertex_distances[thrust::get<0>(t).y] +
                         vertex_distances[thrust::get<0>(t).z])/3.0;
  }
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void DepthSortPolygons(vtkPistonDataObject *id, double *cameravec)
{
  vtkPistonReference *tr = id->GetReference();
  if (tr->type != VTK_POLY_DATA || tr->data == NULL) {
    // Type mismatch, don't bother trying
    return;
  }
  vtk_polydata *pD = (vtk_polydata *)tr->data;

  //
  // we need to compute the distance to the camera for each cell.
  // Perform a dot product of each vertex with the supplied camera vector
  //

  // prepare an array for the distances
  thrust::device_vector<float> distances(pD->points->size());

  // initialize our functor which will compute distance and store in a vector
  float3 cam = make_float3(cameravec[0], cameravec[1], cameravec[2]);
  distance_functor distance(cam);

  // apply distance functor using input and output arrays using zip_iterator
  thrust::for_each(
    thrust::make_zip_iterator(thrust::make_tuple(pD->points->begin(), distances.begin())),
    thrust::make_zip_iterator(thrust::make_tuple(pD->points->end(),   distances.end())),
    distance);

  // to test if it is working, copy the distances into the scalars
  // so we can colour by scalar values
//  thrust::copy(distances.begin(), distances.end(), pD->scalars->begin());
  
  //
  // To compute the average distance for each cell, we must
  // sum/gather 3 distances (one for each vertex) for every cell by
  // looking up the vertex indices from the cell array tuples
  //

  // prepare an array for the distances
  thrust::device_vector<float> cell_distances(pD->nCells);

  celldistance_functor celldist(thrust::raw_pointer_cast(distances.data()));

  thrust::for_each(
    thrust::make_zip_iterator(thrust::make_tuple(pD->cells->begin(), cell_distances.begin())),
    thrust::make_zip_iterator(thrust::make_tuple(pD->cells->end(),   cell_distances.end())),
    celldist);

  //
  // now we want to sort the cells using the average distance
  // we must copy the cell vertex index tuple during the sort
  //
  thrust::sort_by_key(cell_distances.begin(), cell_distances.end(), pD->cells->begin(), 
    thrust::greater<float>());
}

//------------------------------------------------------------------------------
void CudaTransferToGL(vtkPistonDataObject *id, unsigned long dataObjectMTimeCache,
                      vtkPistonScalarsColors *psc,
                      cudaGraphicsResource **vboResources,
                      bool &hasNormals, bool &hasColors, 
                      bool &useindexbuffers)
{
  vtkPistonReference *tr = id->GetReference();
  if (tr->type != VTK_POLY_DATA || tr->data == NULL)
    {
    // Type mismatch, don't bother trying
    return;
    }

  vtk_polydata *pD = (vtk_polydata *)tr->data;

  // Claim access to buffer for cuda
  cudaError_t res;
  res = cudaGraphicsMapResources(4, vboResources, 0);
  if (res != cudaSuccess)
  {
    cerr << "Claim for CUDA failed ... " << cudaGetErrorString(res) << endl;
    return;
  }

  size_t num_bytes;
  float3 *vertexBufferData;
  uint3  *cellsBufferData;
  float  *normalsBufferData;
  // float3 *colorsBufferData;
  float4  *colorbufferdata; 

  res = cudaGraphicsResourceGetMappedPointer
      ((void **)&vertexBufferData, &num_bytes, vboResources[0]);
  if(res != cudaSuccess) {
    cerr << "Get mappedpointer for vertices failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }
  res = cudaGraphicsResourceGetMappedPointer
      ((void **)&normalsBufferData, &num_bytes, vboResources[1]);
  if(res != cudaSuccess) {
    cerr << "Get mappedpointer for normals failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }
  res = cudaGraphicsResourceGetMappedPointer
      ((void **)&colorbufferdata, &num_bytes, vboResources[2]);
  if(res != cudaSuccess)
  {
    cerr << "Get mappedpointer for colors failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }

  res = cudaGraphicsResourceGetMappedPointer
      ((void **)&cellsBufferData, &num_bytes, vboResources[3]);
  if(res != cudaSuccess)
  {
    std::string errormsg = cudaGetErrorString(res);
    cerr << "Get mappedpointer for cell indices failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }

  // Copy on card verts to the shared on card gl buffer
  thrust::copy(pD->points->begin(), pD->points->end(),
               thrust::device_ptr<float3>(vertexBufferData));

  // Copy on card cell indices to the shared on card gl buffer
  if (pD->cells) {
    useindexbuffers = true;
    thrust::copy(pD->cells->begin(), pD->cells->end(),
                 thrust::device_ptr<uint3>(cellsBufferData));
  }

  hasNormals = false;
  if (pD->normals)
    {
    hasNormals = true;

    // Copy on card verts to the shared on card gl buffer
    thrust::copy(pD->normals->begin(), pD->normals->end(),
                 thrust::device_ptr<float>(normalsBufferData));
    }
  hasColors = false;


  if (pD->colors)
  {
//    thrust::fill(pD->colors->begin(), pD->colors->end(), 127);
    thrust::copy(pD->colors->begin(), pD->colors->end(), 
      thrust::device_ptr<float4>(colorbufferdata));
  }
//  else 
/*
  if (pD->scalars)
    {
    double scalarRange[2];
    id->GetScalarsRange(scalarRange);

    hasColors = true;

//    if(id->GetMTime() > dataObjectMTimeCache)
      {
      vtkPiston::minmax_pair<float> result = vtkPiston::find_min_max(
                                              pD->scalars);

      scalarRange[0] = static_cast<double>(result.min_val);
      scalarRange[1] = static_cast<double>(result.max_val);

      // Set parameters to compute scalars colors
      const int numvalues = 256;
      id->SetScalarsRange(scalarRange);
      psc->SetTableRange(scalarRange[0], scalarRange[1]);
      psc->SetNumberOfValues(numvalues);
      }

    std::vector<float> *colors = psc->ComputeScalarsColorsf(VTK_RGB);

    // Copy to GPU
    thrust::device_vector<float> onGPU(colors->begin(), colors->end());
    float *raw_ptr = thrust::raw_pointer_cast(&onGPU[0]);

    // Now run each scalar data through the map to choose a color for it

    // \NOTE: Since GPU most likely going to calculate range using single
    // floating point precision, we may lose precision and hence, we need
    // to check if the range min and max are almost equal
    //TODO: Remove this when piston gives us exactly same values for
    //isocontour.
    float tempRange[2] =
      {
      static_cast<float>(scalarRange[0]),
      static_cast<float>(scalarRange[1])
      };
    if( AlmostEqualRelativeAndAbs(scalarRange[0], scalarRange[1],
                                  numeric_limits<float>::epsilon(),
                                  numeric_limits<float>::epsilon() * 10) )
      {
      tempRange[1] = tempRange[0]+1.0;
      }

    color_map<float> colorMap(raw_ptr, onGPU.size(), VTK_RGB, tempRange[0], tempRange[1]);
    thrust::copy(thrust::make_transform_iterator(pD->scalars->begin(), colorMap),
                 thrust::make_transform_iterator(pD->scalars->end(), colorMap),
                 thrust::device_ptr<float3>(colorsBufferData));
    }
*/
  // Allow GL to access again
  res = cudaGraphicsUnmapResources(4, vboResources, 0);
  if (res != cudaSuccess)
  {
    cerr << "Release from CUDA failed ... " << cudaGetErrorString(res) << endl;
    return;
  }

  return;
}
//------------------------------------------------------------------------------
} //namespace
