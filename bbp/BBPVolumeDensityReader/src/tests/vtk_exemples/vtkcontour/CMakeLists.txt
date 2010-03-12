cmake_minimum_required(VERSION 2.6)
PROJECT(vtkcontour)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(vtkcontour vtkcontour.cxx)
TARGET_LINK_LIBRARIES(vtkcontour vtkHybrid)
