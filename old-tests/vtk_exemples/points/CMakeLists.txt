cmake_minimum_required(VERSION 2.6)
PROJECT(points)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(points points.cxx)
TARGET_LINK_LIBRARIES(points vtkHybrid)
