cmake_minimum_required(VERSION 2.6)
PROJECT(vtkImageAlgorithmFilter)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
 
ADD_EXECUTABLE(vtkImageAlgorithmFilter main.cxx vtkImageAlgorithmFilter.cxx)
TARGET_LINK_LIBRARIES(vtkImageAlgorithmFilter vtkVolumeRendering vtkHybrid)
