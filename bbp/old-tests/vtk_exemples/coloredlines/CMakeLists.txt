cmake_minimum_required(VERSION 2.6)
PROJECT(coloredlines)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(coloredlines coloredlines.cxx)
TARGET_LINK_LIBRARIES(coloredlines vtkHybrid)
