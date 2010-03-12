cmake_minimum_required(VERSION 2.6)
PROJECT(pickingpixel)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(pickingpixel pickingpixel.cxx)
TARGET_LINK_LIBRARIES(pickingpixel vtkHybrid)
