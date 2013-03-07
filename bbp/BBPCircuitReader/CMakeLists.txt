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
include(vtkModuleAPI)
include(vtkModuleMacros)
include(ParaViewMacros)

find_package(Qt4 )

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
# pv-zoltan
#--------------------------------------------------
set(ZOLTAN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../pv-zoltan")
if (EXISTS "${ZOLTAN_DIR}" AND IS_DIRECTORY "${ZOLTAN_DIR}")
  set(ZOLTAN_SUPPORT 1)
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../pv-zoltan
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../pv-zoltan/trilinos/packages/zoltan/src/include
    ${CMAKE_CURRENT_BINARY_DIR}/../../../pv-zoltan/trilinos/packages/zoltan/src
  )
  if(WIN32)
    include_directories(
      ${CMAKE_CURRENT_SOURCE_DIR}/../../../pv-zoltan/trilinos/cmake/tribits/common_tools/win_interface/include
   )
  endif(WIN32)
  ADD_DEFINITIONS(-DZOLTAN_SUPPORT)
  set(ZOLTAN_LIBS pv_zoltan)
  message("Zoltan support enabled")
endif (EXISTS "${ZOLTAN_DIR}" AND IS_DIRECTORY "${ZOLTAN_DIR}")


#--------------------------------------------------
# Includes
#--------------------------------------------------
include_directories(
  ${PLUGIN_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/Piston
  ${BBP_SDK_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIR_CPP}
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
  REQUIRED_PLUGINS pv_zoltan

  SERVER_MANAGER_XML     
    ${XML_PATH}/vtkTemporalDifferenceFilter.xml
    ${XML_PATH}/vtkExponentialDecayFilter.xml
    ${XML_PATH}/vtkNeuronAlphaFunction.xml
    #
    ${XML_PATH}/vtkCircuitReader.xml
    ${XML_PATH}/vtkMorphologyReader.xml
    ${XML_PATH}/vtkMeshBinReader.xml
    ${XML_PATH}/vtkH5MCellReader.xml

  SERVER_MANAGER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalDifferenceFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkExponentialDecayFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkNeuronAlphaFunction.cxx
    #
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCircuitReader.cxx
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
    pv_zoltan
    ${BBPSDK_LIBRARIES}
    hdf5_cpp
    hdf5
    ${CUDA_LIBRARIES}
#    ${BBPSDK_LIBRARIES}
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

message("BBPSDK LIBS ARE 
    ${BBPSDK_LIBRARIES}
")

message("BOOST LIBS ARE 
    ${Boost_LIBRARIES}
")

message("HDF5 LIBS ARE 
    ${HDF5_C_LIBRARY} 
    ${HDF5_CXX_LIBRARY} 
")


add_executable(test_bbplink test_bbplink.cpp)
target_link_libraries(test_bbplink 
    ${PLUGIN_NAME}
#    ${BBPSDK_LIBRARIES}
#    ${CUDA_LIBRARIES}
#    ${ZOLTAN_LIBS}
#    ${HDF5_C_LIBRARY} 
#    ${HDF5_CXX_LIBRARY} 
#    ${Boost_LIBRARIES}


)

