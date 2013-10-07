cmake_minimum_required(VERSION 2.6)
PROJECT(coloredlines)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(simple_vtkImageData simple_vtkImageData.cxx)
TARGET_LINK_LIBRARIES(simple_vtkImageData vtkHybrid)
