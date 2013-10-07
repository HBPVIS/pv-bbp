cmake_minimum_required(VERSION 2.6)
PROJECT(polylines)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(polylines polylines.cxx)
TARGET_LINK_LIBRARIES(polylines vtkHybrid)
