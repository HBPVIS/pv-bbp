#
# 
#       Ecole Polytechnique Federale de Lausanne
#       Brain Mind Institute,
#       Blue Brain Project
#       (c) 2005-2010. All rights reserved.
# 
#       Authors: Sebastien Lasserre
#

# create a paraview plugin containing server manager xml and the server
# manager classes to build
# this plugin can be loaded on the server side

project(BBPCircuitReader)
cmake_minimum_required(VERSION 2.8)

# Policy issue with cmake 
# FIXME : there are still some warnings during the cmake process
if(COMMAND cmake_policy)
     cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# plugin path and CMake module path
get_filename_component(PLUGIN_PATH ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)

# Locate ParaView build and then import CMake configuration, 
# macros etc. from it.
# FIXME : currently the user need to run cmake with the VTK_DIR definition :
# ccmake -DVTK_DIR=/home/lasserre/local/ParaView3/build/VTK/
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
#INCLUDE(${PARAVIEW_VTK_DIR}/VTKConfig.cmake)

find_package(Qt4 )

# --------------------------------------------------
# Optional Trilinos support 
# --------------------------------------------------
SET(TRILINOS_ENABLED 0)
IF (PARAVIEW_USE_MPI)
  OPTION(USE_TRILINOS "Use Trilinos library" OFF)
  IF(USE_TRILINOS)

    OPTION(USE_SYSTEM_TRILINOS "Use System Trilinos library" OFF)
    IF(USE_SYSTEM_TRILINOS)
      FIND_PACKAGE(Trilinos)
    ELSE(USE_SYSTEM_TRILINOS)
      SET(TPL_ENABLE_MPI ON CACHE BOOL "Do not change")
      SET(Trilinos_ENABLE_TESTS OFF CACHE BOOL "Do not change")
      SET(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF CACHE BOOL "Do not change")
      SET(Trilinos_ENABLE_Zoltan ON CACHE BOOL "Do not change")
      SET(Trilinos_ENABLE_Fortran OFF CACHE BOOL "Do not change")
      IF (WIN32)
        ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
        INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/trilinos/cmake/tribits/common_tools/win_interface/include")
      ENDIF(WIN32)
      ADD_SUBDIRECTORY(trilinos)
      SET(Trilinos_FOUND 1)
      SET(Trilinos_DEPS_HTML_OUTPUT_FILE "" CACHE STRING "Leave empty" FORCE)
    ENDIF(USE_SYSTEM_TRILINOS)

    IF(Trilinos_FOUND)
      INCLUDE_DIRECTORIES("${Trilinos_INCLUDE_DIRS}")
      LINK_DIRECTORIES("${Trilinos_LIBRARY_DIRS}")
      # MESSAGE("Trilinos Include Dirs ${Trilinos_INCLUDE_DIRS} and Trilinos Link Dirs ${Trilinos_LIBRARY_DIRS} ")
      #
      # Add source which depends on trilinos to our plugin
      #
      SET(TRILINOS_CXX ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.cxx)
      SET(TRILINOS_XML ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.xml)
      SET(TRILINOS_LIBS zoltan)
      SET(TRILINOS_ENABLED 1)
      ADD_DEFINITIONS(-DPLUGIN_USE_TRILINOS)
    ENDIF(Trilinos_FOUND)

  ENDIF(USE_TRILINOS)
ELSE (PARAVIEW_USE_MPI)
  SET(USE_TRILINOS OFF CACHE BOOL "Use Trilinos library" FORCE)
ENDIF (PARAVIEW_USE_MPI)

# --------------------------------------------------
# Plugin config
# --------------------------------------------------
# PLUGIN PATHS
set(XML_PATH   "${CMAKE_CURRENT_SOURCE_DIR}")
set(SRC_PATH   "${CMAKE_CURRENT_SOURCE_DIR}")

# INCLUDES
include_directories(
  ${PLUGIN_PATH}
  ${BBP_SDK_INCLUDE_DIRS}
)


#--------------------------------------------------
# CONDITIONAL : Setup Qt/GUI panel sources/wrapping
#--------------------------------------------------
IF(PARAVIEW_BUILD_QT_GUI)

  QT4_WRAP_CPP(
    BBP_MOC_SRCS 
      ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanel.h
  )

  ADD_PARAVIEW_OBJECT_PANEL(
      BBP_IFACE 
      BBP_IFACE_SRCS 
    CLASS_NAME 
      pqCircuitReaderPanel
    XML_NAME 
      CircuitReader 
    XML_GROUP 
      sources
  )

  QT4_WRAP_UI(
    BBP_UI_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanel.ui
  )

ENDIF(PARAVIEW_BUILD_QT_GUI)

#--------------------------------------------------
# Define plugin sources
#--------------------------------------------------

SET(PLUGIN_NAME pv_BBP)
ADD_PARAVIEW_PLUGIN(
  ${PLUGIN_NAME}
  "1.0" 

  SERVER_MANAGER_XML     
    ${XML_PATH}/vtkCircuitReader.xml
    ${XML_PATH}/vtkTemporalDifferenceFilter.xml
    ${TRILINOS_XML}
    ${XML_PATH}/vtkPartitionOutline.xml
    ${XML_PATH}/vtkDepthSortRepresentation.xml
    ${XML_PATH}/vtkMorphologyReader.xml
    ${XML_PATH}/vtkMeshBinReader.xml
    ${XML_PATH}/vtkH5MCellReader.xml

  SERVER_MANAGER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCircuitReader.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalDifferenceFilter.cxx
    ${TRILINOS_CXX}
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkBoundsExtentTranslator.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkPartitionOutline.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkDepthSortRepresentation.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkDepthSortDefaultPainter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkDepthSortPainter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTwoScalarsToColorsPainter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMorphologyReader.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMeshBinReader.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5MCellReader.cxx

  SERVER_SOURCES

  GUI_INTERFACES 
    ${BBP_IFACE}

  GUI_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanel.cxx
    ${BBP_IFACE_SRCS}
    ${BBP_MOC_SRCS}
    ${BBP_UI_SOURCES}

  GUI_RESOURCE_FILES
    ${XML_PATH}/pv_BBP_Readers.xml
    ${XML_PATH}/pv_BBP_Sources.xml
    ${XML_PATH}/pv_BBP_Filters.xml
)

#--------------------------------------------------
# setup plugin linking
#--------------------------------------------------
TARGET_LINK_LIBRARIES(${PLUGIN_NAME} 
        vtkRenderingFreeTypeFontConfig
    ${BBP_SDK_LIBRARIES}
    ${TRILINOS_LIBS}
    ${HDF5_C_LIBRARY} 
    ${HDF5_CXX_LIBRARY} 
    ${Boost_SERIALIZATION_LIBRARY} 
    ${Boost_THREAD_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_REGEX_LIBRARY}
    ${Boost_SERIALIZATION_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    ${Boost_DATE_TIME_LIBRARY}
)


