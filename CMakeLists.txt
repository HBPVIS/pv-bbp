cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(BrainView)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake
  ${CMAKE_SOURCE_DIR}/CMake/common)

include(GitExternal)

set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "0")
set(VERSION_ABI 1)

set(BRAINVIEW_DESCRIPTION "Blue Brain Plugin for ParaView")
set(BRAINVIEW_MAINTAINER "CSCS <help@cscs.ch>")
set(BRAINVIEW_LICENSE BSD)

set(BUILD_TESTING OFF CACHE BOOL "Tests")
set(BUILD_EXAMPLES OFF CACHE BOOL "Examples")
set(PARAVIEW_QT_VERSION 5 CACHE STRING "Qt Version")
set(PARAVIEW_USE_MPI ON CACHE BOOL "Use MPI")

set(COMMON_PROJECT_DOMAIN ch.epfl.bluebrain)
include(Common)

common_package(Boost REQUIRED COMPONENTS filesystem program_options
  serialization system thread unit_test_framework)
common_package(HDF5 REQUIRED)
common_package(BBPSDK REQUIRED)
common_package(ParaView REQUIRED)
common_package(pvZoltan)

include(${PARAVIEW_USE_FILE})
include(vtkModuleAPI)
include(vtkModuleMacros)
include(ParaViewMacros)

# --------------------------------------------------------------
# setup paths to 'default' data
# --------------------------------------------------------------
site_name(MACHINE)
if (MACHINE MATCHES "carona")
  set(PV_BBP_DEFAULT_DATA "/Users/biddisco/data/bbp/egpgv/centralV.cfg")
  set(PV_BBP_DEFAULT_TARGET "5K")
elseif (MACHINE MATCHES "bbpviz")
  set (PV_BBP_DEFAULT_DATA "/gpfs/bbp.cscs.ch/project/proj3/resources/circuits/summit2015demo/BlueConfig")
  set(PV_BBP_DEFAULT_TARGET "demo")
elseif (MACHINE MATCHES "crusca")
  set (PV_BBP_DEFAULT_DATA "/gpfs/bbp.cscs.ch/project/proj3/resources/circuits/summit2015demo/BlueConfig")
  set(PV_BBP_DEFAULT_TARGET "1K")
else ()
  set (PV_BBP_DEFAULT_DATA "/gpfs/bbp.cscs.ch/project/proj3/resources/circuits/summit2015demo/BlueConfig")
  set(PV_BBP_DEFAULT_TARGET "demo")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkCircuitReaderBase.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/vtkCircuitReaderBase.xml
  @ONLY
  )

# --------------------------------------------------------------
# Build main plugin
# --------------------------------------------------------------

# --------------------------------------------------
# Plugin config
# --------------------------------------------------
# PLUGIN PATHS
set(XML_PATH   "${CMAKE_CURRENT_SOURCE_DIR}")
set(SRC_PATH   "${CMAKE_CURRENT_SOURCE_DIR}")

#--------------------------------------------------
# pv-zoltan : Required for depth sort
#--------------------------------------------------
if(PVZOLTAN_FOUND)
  set(ZOLTAN_SUPPORT 1)
  ADD_DEFINITIONS(-DZOLTAN_SUPPORT -DPV_BBP_USE_ZOLTAN)
  set(ZOLTAN_LIBS pv_zoltan)
  include_directories(${COMMON_SOURCE_DIR}/pvZoltan ${Trilinos_INCLUDE_DIRS})
  set(PVZOLTAN_TARGET pvZoltan)
endif()

#--------------------------------------------------
# pv-temporal : Required for spike display
#--------------------------------------------------
set(TEMPORAL_SUPPORT 1)
ADD_DEFINITIONS(-DTEMPORAL_SUPPORT -DPV_BBP_USE_TEMPORAL)
set(TEMPORAL_LIBS pv_temporal)
include_directories(${COMMON_SOURCE_DIR}/pv-temporal)

#--------------------------------------------------
# Includes
#--------------------------------------------------
get_filename_component(PLUGIN_PATH ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)
include_directories(
  ${PLUGIN_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/Piston
  )

#--------------------------------------------------
# CONDITIONAL : Setup Qt/GUI panel sources/wrapping
#--------------------------------------------------
IF(PARAVIEW_BUILD_QT_GUI)
  if (PARAVIEW_QT_VERSION VERSION_GREATER "4")
    set (Qt5_FIND_COMPONENTS Widgets)
    include (ParaViewQt5)

    QT5_WRAP_CPP(
      BBP_MOC_SRCS
      ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelBase.h
      ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelMesh.h
      ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelSoma.h
      )

    QT5_WRAP_UI(
      BBP_UI_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelBase.ui
      )

  else ()
    include (${QT_USE_FILE})

    QT4_WRAP_CPP(
      BBP_MOC_SRCS
      ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelBase.h
      ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelMesh.h
      ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelSoma.h
      )

    QT4_WRAP_UI(
      BBP_UI_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelBase.ui
      )
  endif ()

  ADD_PARAVIEW_OBJECT_PANEL(
    BBP_IFACE
    BBP_IFACE_SRCS
    CLASS_NAME
    pqCircuitReaderPanelMesh
    XML_NAME
    CircuitReaderMesh
    XML_GROUP
    sources
    )

  ADD_PARAVIEW_OBJECT_PANEL(
    BBP_IFACE
    BBP_IFACE_SRCS
    CLASS_NAME
    pqCircuitReaderPanelSoma
    XML_NAME
    CircuitReaderSoma
    XML_GROUP
    sources
    )
ENDIF(PARAVIEW_BUILD_QT_GUI)

#--------------------------------------------------
# Define plugin sources
#--------------------------------------------------

SET(PLUGIN_NAME pv_BBP)
ADD_PARAVIEW_PLUGIN(
  ${PLUGIN_NAME}
  "1.0"
  REQUIRED_PLUGINS ${PVZOLTAN_TARGET} pv_temporal

  SERVER_MANAGER_XML
  ${XML_PATH}/vtkNeuronAlphaFunction.xml
  #
  ${CMAKE_CURRENT_BINARY_DIR}/vtkCircuitReaderBase.xml
  ${XML_PATH}/vtkCircuitReaderMesh.xml
  ${XML_PATH}/vtkCircuitReaderSoma.xml
  ${XML_PATH}/vtkNeuronSpikeFilter.xml
  ${XML_PATH}/vtkNeuronSpikeTableSource.xml

  SERVER_MANAGER_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkNeuronAlphaFunction.cxx
  #
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkCircuitReaderBase.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkCircuitReaderMesh.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkCircuitReaderSoma.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkNeuronSpikeFilter.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkNeuronSpikeTableSource.cxx

  SERVER_SOURCES

  GUI_INTERFACES
  ${BBP_IFACE}

  GUI_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelBase.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelMesh.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/pqCircuitReaderPanelSoma.cxx
  ${BBP_IFACE_SRCS}
  ${BBP_MOC_SRCS}
  ${BBP_UI_SOURCES}
  )

if(ADD_CXX)
  #  message("Adding C++0x flag for Depth Sort")
  #  set_target_properties(pv_BBP PROPERTIES COMPILE_FLAGS "-std=c++0x -I/usr/include/c++/4.2.1/tr1")
endif()

#--------------------------------------------------
# setup plugin linking
#--------------------------------------------------

TARGET_LINK_LIBRARIES(${PLUGIN_NAME}
  PUBLIC
  ${PVZOLTAN_TARGET}
  pv_temporal
  vtksys
  BBPSDK Brion Lunchbox pthread
  hdf5_cpp
  hdf5
  ${CUDA_LIBRARIES}
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

#--------------------------------------------------
# A simple test of link against plugin
#--------------------------------------------------
add_executable(test_bbplink test_bbplink.cpp)
target_link_libraries(test_bbplink
  ${PLUGIN_NAME}
  )

return()

#--------------------------------------------------
# Install
#--------------------------------------------------
SET(INSTALL_PATH
  "${CMAKE_INSTALL_PREFIX}/lib/paraview-${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}"
  )

INSTALL(
  TARGETS
  ${PLUGIN_NAME}
  DESTINATION ${INSTALL_PATH}
  )
