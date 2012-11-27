cmake_minimum_required(VERSION 2.6)
PROJECT(picking)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(picking picking.cxx)
TARGET_LINK_LIBRARIES(picking vtkHybrid)
