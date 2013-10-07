cmake_minimum_required(VERSION 2.6)
PROJECT(volrendering)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(volrendering volrendering.cxx)
TARGET_LINK_LIBRARIES(volrendering vtkHybrid vtkVolumeRendering)
