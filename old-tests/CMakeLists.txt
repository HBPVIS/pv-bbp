cmake_minimum_required(VERSION 2.6)
PROJECT(VisualizeMorphology)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# BBP-SDK
set(CMAKE_MODULE_PATH   "${CMAKE_SOURCE_DIR}/Modules")
FIND_PACKAGE(BBP-SDK REQUIRED)
 
# INCLUDES
include_directories(${BBP-SDK_INCLUDE_DIRS})

ADD_EXECUTABLE(VisualizeMorphology VisualizeMorphology.cxx vtkMorphologyReader.cxx)
TARGET_LINK_LIBRARIES(VisualizeMorphology vtkHybrid ${BBP-SDK_LIB})
