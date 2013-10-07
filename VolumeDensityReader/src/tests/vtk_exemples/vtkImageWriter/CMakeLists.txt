cmake_minimum_required(VERSION 2.6)
PROJECT(vtkimagewriter)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(vtkimagewriter vtkimagewriter.cxx)
TARGET_LINK_LIBRARIES(vtkimagewriter vtkHybrid)
