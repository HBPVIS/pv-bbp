cmake_minimum_required(VERSION 2.6)
PROJECT(TubesWithVaryingRadiusAndColors)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(TubesWithVaryingRadiusAndColors TubesWithVaryingRadiusAndColors.cxx)
TARGET_LINK_LIBRARIES(TubesWithVaryingRadiusAndColors vtkHybrid)
