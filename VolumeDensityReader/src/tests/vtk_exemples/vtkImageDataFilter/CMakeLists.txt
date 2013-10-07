cmake_minimum_required(VERSION 2.6)
PROJECT(vtkImageProcess)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
 
ADD_EXECUTABLE(vtkImageProcess main.cxx vtkImageProcess.cxx)
TARGET_LINK_LIBRARIES(vtkImageProcess vtkVolumeRendering vtkHybrid)
