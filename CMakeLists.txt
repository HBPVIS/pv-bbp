project(ParaBBP)

cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
FIND_PACKAGE(BBP-SDK REQUIRED)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bbp/BBPVolumeDensityReader/proj/cmake)
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