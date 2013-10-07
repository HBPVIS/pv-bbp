cmake_minimum_required(VERSION 2.6)
PROJECT(simple_vtkImageData_viz)
 

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# BBP-SDK
set(CMAKE_MODULE_PATH   "${CMAKE_SOURCE_DIR}/Modules")
FIND_PACKAGE(BBP-SDK REQUIRED)

# INCLUDES
include_directories(${BBP-SDK_INCLUDE_DIRS})


ADD_EXECUTABLE(simple_vtkImageData_viz simple_vtkImageData_viz.cxx)
TARGET_LINK_LIBRARIES(simple_vtkImageData_viz vtkHybrid vtkVolumeRendering ${BBP-SDK_LIB})
