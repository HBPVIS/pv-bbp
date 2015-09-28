project(pv-BBP)
cmake_minimum_required(VERSION 2.8)

# --------------------------------------------------------------
# CMake setup, modules, policy, paths
# --------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
  cmake_policy(SET CMP0017 NEW)
  cmake_policy(SET CMP0023 OLD)
endif(COMMAND cmake_policy)
      
# CMake 2.6+ specific
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# --------------------------------------------------------------
# Plugin path
# --------------------------------------------------------------
get_filename_component(PLUGIN_PATH ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)

# --------------------------------------------------------------
# Find BBP-SDK
# --------------------------------------------------------------
find_package(BBPSDK REQUIRED NO_MODULE)
include_directories(${BBPSDK_INCLUDE_DIRS})
add_definitions(-DCXX_UNORDERED_MAP_SUPPORTED)

# --------------------------------------------------------------
# Find ParaView
# --------------------------------------------------------------
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
# recent VTK/ParaView build uses special module files
include(vtkModuleAPI)
include(vtkModuleMacros)
include(ParaViewMacros)

# --------------------------------------------------------------
# Find Boost
# --------------------------------------------------------------
OPTION(Boost_USE_FIND_NO_MODULE "If cmake installed boost, use NO_MODULE in find package" OFF)

if(Boost_USE_FIND_NO_MODULE)
  find_package(Boost NO_MODULE COMPONENTS program_options;filesystem;system;regex;thread;date_time;serialization)
else(Boost_USE_FIND_NO_MODULE)
  find_package(Boost 1.41.0 REQUIRED filesystem program_options serialization system thread unit_test_framework)
endif(Boost_USE_FIND_NO_MODULE)

include_directories(
  ${Boost_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
)
link_directories(${Boost_LIBRARY_DIRS})

# --------------------------------------------------------------
# Find HDF5
# --------------------------------------------------------------
find_package(HDF5 REQUIRED NO_MODULE PATHS "$ENV{ProgramFiles}/HDF5")
include_directories(${HDF5_INCLUDE_DIR}
                    ${HDF5_INCLUDE_DIR_CPP})

# --------------------------------------------------------------
# set c++0x flag for bbp #includes
# --------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(ADD_CXX 1)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(ADD_CXX 1)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
endif()

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
  set(PV_BBP_DEFAULT_TARGET "1K")
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
# HDF5
#--------------------------------------------------
FIND_PACKAGE(HDF5 REQUIRED NO_MODULE PATHS "$ENV{ProgramFiles}/HDF5")
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR} ${HDF5_INCLUDE_DIR_CPP})

#--------------------------------------------------
# pv-zoltan : Required for depth sort
#--------------------------------------------------
find_package(pv_zoltan NO_MODULE REQUIRED)
if(NOT pv_zoltan_FOUND)
  message(FATAL_ERROR "pv-bbp requires pv_zoltan")
endif()
include_directories(${pv_zoltan_INCLUDE_DIR})

set(ZOLTAN_SUPPORT 1)
ADD_DEFINITIONS(-DZOLTAN_SUPPORT -DPV_BBP_USE_ZOLTAN)
set(ZOLTAN_LIBS pv_zoltan)
message("Zoltan support enabled")

#--------------------------------------------------
# pv-temporal : Required for spike display
#--------------------------------------------------
find_package(pv_temporal NO_MODULE REQUIRED)
if(NOT pv_temporal_FOUND)
  message(FATAL_ERROR "pv-bbp requires pv_temporal")
endif()
include_directories(${pv_temporal_INCLUDE_DIR})

set(TEMPORAL_SUPPORT 1)
ADD_DEFINITIONS(-DTEMPORAL_SUPPORT -DPV_BBP_USE_TEMPORAL)
set(TEMPORAL_LIBS pv_temporal)
message("Temporal support enabled")

#--------------------------------------------------
# Includes
#--------------------------------------------------
include_directories(
  ${PLUGIN_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/Piston
  ${HDF5_INCLUDE_DIR_CPP}
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
  REQUIRED_PLUGINS pv_zoltan pv_temporal

  SERVER_MANAGER_XML     
    ${XML_PATH}/vtkNeuronAlphaFunction.xml
    #
    ${CMAKE_CURRENT_BINARY_DIR}/vtkCircuitReaderBase.xml
    ${XML_PATH}/vtkCircuitReaderMesh.xml
    ${XML_PATH}/vtkCircuitReaderSoma.xml
    ${XML_PATH}/vtkNeuronSpikeFilter.xml

  SERVER_MANAGER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkNeuronAlphaFunction.cxx
    #
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCircuitReaderBase.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCircuitReaderMesh.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCircuitReaderSoma.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkNeuronSpikeFilter.cxx

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
    pv_zoltan
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
