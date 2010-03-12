cmake_minimum_required(VERSION 2.6)
PROJECT(vtkimageio)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(vtkimageio vtkimageio.cxx)
TARGET_LINK_LIBRARIES(vtkimageio vtkHybrid)
