cmake_minimum_required(VERSION 2.6)
PROJECT(VisualizeVolumeDensity)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# BBP-SDK
set(CMAKE_MODULE_PATH   "${CMAKE_SOURCE_DIR}/Modules")
FIND_PACKAGE(BBP-SDK REQUIRED)
 
# INCLUDES
include_directories(${BBP-SDK_INCLUDE_DIRS})

ADD_EXECUTABLE(VisualizeVolumeDensity VisualizeVolumeDensity.cxx vtkVolumeDensityReader.cxx)
TARGET_LINK_LIBRARIES(VisualizeVolumeDensity vtkVolumeRendering vtkHybrid)
TARGET_LINK_LIBRARIES(VisualizeVolumeDensity vtkVolumeRendering vtkHybrid ${BBP-SDK_LIB})
