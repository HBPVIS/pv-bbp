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

project(BBPMeshBinReader)
cmake_minimum_required(VERSION 2.8)

# Policy issue with cmake 
# FIXME : there are still some warnings during the cmake process
if(COMMAND cmake_policy)
     cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# plugin path and CMake module path
get_filename_component(PLUGIN_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../.. ABSOLUTE)

# Locate ParaView build and then import CMake configuration, 
# macros etc. from it.
# FIXME : currently the user need to run cmake with the VTK_DIR definition :
# ccmake -DVTK_DIR=/home/lasserre/local/ParaView3/build/VTK/
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
INCLUDE(${PARAVIEW_VTK_DIR}/VTKConfig.cmake)

# PLUGIN PATHS
set(XML_PATH   "${PLUGIN_PATH}/xml")
set(SRC_PATH   "${PLUGIN_PATH}/src")

# INCLUDES
include_directories(
    ${PLUGIN_PATH}/src
    ${BBP_SDK_INCLUDE_DIRS})

#--------------------------------------------------
# Define Plugin 
#--------------------------------------------------
ADD_PARAVIEW_PLUGIN(
  BBPMeshBinReader 
  "1.0"
  SERVER_MANAGER_XML     ${XML_PATH}/MeshBinReader.xml 
  SERVER_MANAGER_SOURCES ${SRC_PATH}/vtkMeshBinReader.cxx
  GUI_RESOURCE_FILES     ${XML_PATH}/GUIMeshBinReader.xml
  )

# LINK
TARGET_LINK_LIBRARIES(BBPMeshBinReader 
    ${BBP_SDK_LIBRARIES}
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

