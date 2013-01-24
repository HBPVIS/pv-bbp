project(pv-BBP)

cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# --------------------------------------------------------------
# Find BBP-SDK
# --------------------------------------------------------------
SET(BBP_SDK_USE_STATIC_LIBS ON)
find_package(BBP-SDK 0.14 REQUIRED common io model)
include_directories(${BBP_SDK_INCLUDE_DIRS}  ${BBP_SDK_PATH}/include)

message("BBP_SDK_INCLUDE_DIRS are ${BBP_SDK_INCLUDE_DIRS}")


# --------------------------------------------------------------
# Find Boost
# --------------------------------------------------------------
IF (MSVC)
	message("Setting MSVC boost options")
	set(Boost_USE_STATIC_LIBS        ON)
	set(Boost_USE_STATIC             ON)
	set(Boost_USE_MULTITHREADED      ON)
	set(Boost_USE_STATIC_RUNTIME    OFF)
ENDIF (MSVC)

#find_package(Boost 1.49 REQUIRED
#  thread
#  filesystem
#  system
#  serialization
#  program_options
#  date_time
#  regex
#  unit_test_framework
#)
#link_directories(${Boost_LIBRARY_DIRS})
#include_directories(AFTER SYSTEM ${Boost_INCLUDE_DIRS})

# --------------------------------------------------------------
# Find HDF5
# --------------------------------------------------------------
#find_package(HDF5 COMPONENTS C CXX REQUIRED NO_MODULE PATHS "$ENV{ProgramFiles}/HDF5")
#include_directories(${HDF5_INCLUDE_DIR}
#                    ${HDF5_INCLUDE_DIR}/cpp)

#--------------------------------------------------
# If (paraview) using system H5FDdsm
#--------------------------------------------------
IF (NOT TARGET "H5FDdsm")
#  FIND_PACKAGE(H5FDdsm REQUIRED NO_MODULE PATHS "$ENV{ProgramFiles}/H5FDdsm")
#  INCLUDE_DIRECTORIES(${H5FDdsm_INCLUDE_DIR} ${H5FD_DSM_INCLUDE_DIR})
ENDIF (NOT TARGET "H5FDdsm")

# CMake 2.6+ specific
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPCircuitReader)
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/proj/cmake)
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPMeshBinReader/proj/cmake)
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPMorphologyReader/proj/cmake)


#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPMeshBinReader/src/tests/viz_mesh/CMakeLists.txt

#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPMorphologyReader/src/tests/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPMorphologyReader/src/tests/vtk_exemples/coloredlines/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPMorphologyReader/src/tests/vtk_exemples/picking/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPMorphologyReader/src/tests/vtk_exemples/points/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPMorphologyReader/src/tests/vtk_exemples/polylines/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPMorphologyReader/src/tests/vtk_exemples/tubes/CMakeLists.txt


#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/src/tests/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/src/tests/vtk_exemples/pickingpixel/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/src/tests/vtk_exemples/simple_vtkImageData/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/src/tests/vtk_exemples/simple_vtkImageData_viz/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/src/tests/vtk_exemples/volrendering/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/src/tests/vtk_exemples/vtkcontour/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/src/tests/vtk_exemples/vtkImageAlgorithmFilter/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/src/tests/vtk_exemples/vtkImageDataFilter/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/src/tests/vtk_exemples/vtkImageIO/CMakeLists.txt
#add_subdir(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/src/tests/vtk_exemples/vtkImageWriter/CMakeLists.txt
