# Ecole Polytechnique Federale de Lausanne
# Brain Mind Institute,
# Blue Brain Project
# (c) 2006-2011. All rights reserved.
#
# Author: Thomas Traenkler, Barthelemy von Haller
#
#
# - Try to find BBP-SDK
# If BBP_SDK_ROOT is defined this will be used to locate the installation.
# If BBP_SDK_USE_STATIC_LIBS is defined and set to ON it will link against the static library of SDK (hope you compiled it)
#
# Once done this script will define
#  BBP_SDK_FOUND            System has BBP-SDK
#  BBP_SDK_INCLUDE_DIRS     The BBP-SDK include directories
#  BBP_SDK_LIBRARIES        The libraries needed to use BBP-SDK
#  BBP_SDK_DEFINITIONS      Compiler switches required for using BBP-SDK
#  BBP_SDK_DEB_DEPENDENCIES Dependencies for downstream Debian packages
#  (additional components-based variables, see below)
#
# Components
# The components can be : common, io, model, analysis, corba
#
# Be careful to ask for all the components you need plus dependencies. This means
# that asking for model without io and common might lead to problems. 
# TODO: Handle these dependencies.
#
# For each component you specify in find_package(), the following (UPPER-CASE)
# variables are set.  You can use these variables if you would like to pick and
# choose components for your targets instead of just using BBP_SDK_LIBRARIES.
#
#  BBP_SDK_${COMPONENT}_FOUND       True IF the SDK library "component" was found.
#  BBP_SDK_${COMPONENT}_LIBRARY     Contains the libraries for the specified SDK "component" 

cmake_minimum_required(VERSION 2.6)

if(POLICY CMP0011)
  cmake_policy(SET CMP0011 NEW)
endif()

include(FindPackageHandleStandardArgs)

# Versioning
if(BBP-SDK_FIND_VERSION VERSION_LESS 0.11)
  set(BBP-SDK_FIND_SOVERSION 0)
elseif(BBP-SDK_FIND_VERSION VERSION_LESS 0.12)
  set(BBP-SDK_FIND_SOVERSION 1)
elseif(BBP-SDK_FIND_VERSION VERSION_LESS 0.13)
  set(BBP-SDK_FIND_SOVERSION 2)
elseif(BBP-SDK_FIND_VERSION VERSION_LESS 0.14)
  set(BBP-SDK_FIND_SOVERSION 3)
else()
  # TODO: get this from version.h or any other file!
  set(BBP-SDK_SOVERSION 4)
endif()

# Support preference of static libs by adjusting CMAKE_FIND_LIBRARY_SUFFIXES
if(BBP_SDK_USE_STATIC_LIBS)
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a )
  endif()
endif()

# List of the valid SDK components
set(BBP-SDK_VALID_COMPONENTS common io model analysis corba)
if(NOT BBP-SDK_FIND_COMPONENTS)
  set(BBP-SDK_FIND_COMPONENTS common io model)
endif()

if(BBP-SDK_FIND_REQUIRED)
  set(BBP_SDK_VERSION_OUTPUT_TYPE FATAL_ERROR)
else()
  set(BBP_SDK_VERSION_OUTPUT_TYPE STATUS)
endif()

# locate BBP-SDK installation root directory
if(NOT BBP_SDK_ROOT)
  # if defined use environment variable BBP_SDK_ROOT
  set(BBP_SDK_ROOT "$ENV{BBP_SDK_ROOT}")
endif()
# else heuristic search in standard places
if(BBP_SDK_ROOT STREQUAL "")
  set(BBP_SDK_ROOT)
  set(BBP_SDK_SEARCH_PATHS
    ${CMAKE_SOURCE_DIR}/../BBP-SDK/${CMAKE_BUILD_TYPE}/install
    ${CMAKE_SOURCE_DIR}/../../BBP-SDK/${CMAKE_BUILD_TYPE}/install
    ${CMAKE_SOURCE_DIR}/../BBP-SDK/Debug/install
    ${CMAKE_SOURCE_DIR}/../../BBP-SDK/Debug/install
    ${CMAKE_SOURCE_DIR}/../BBP-SDK/Release/install
    ${CMAKE_SOURCE_DIR}/../../BBP-SDK/Release/install
    /usr
    /usr/local
    /opt/BBP-SDK
    /opt/BBP-SDK/latest
    "C:/Program Files/BBP-SDK"
    /opt
    /opt/local
    "D:/cmakebuild/buildyard/install/lib"
  )
  find_path(BBP_SDK_ROOT include/BBP${BBP-SDK_FIND_SOVERSION}/common.h
    NO_DEFAULT_PATH PATHS ${BBP_SDK_SEARCH_PATHS})

  # search again in latest version
  if(NOT BBP_SDK_ROOT)
    find_path(BBP_SDK_ROOT include/BBP/bbp.h
      NO_DEFAULT_PATH PATHS ${BBP_SDK_SEARCH_PATHS})
  endif()
endif()

if(NOT BBP_SDK_ROOT)
    set(BBP_SDK_FAIL TRUE)
    message( ${BBP_SDK_VERSION_OUTPUT_TYPE}
        "ERROR: Can't find BBP-SDK header file bbp.h. Please provide a valid BBP_SDK_ROOT.")
endif()

# find BBP-SDK include directory
unset(BBP_SDK_INCLUDE_DIR CACHE)
find_path(BBP_SDK_INCLUDE_DIR BBP${BBP-SDK_FIND_SOVERSION}/Common/config.h
  HINTS ${BBP_SDK_ROOT}/include)
if(NOT BBP_SDK_INCLUDE_DIR)
  # search again in latest version
  find_path(BBP_SDK_INCLUDE_DIR BBP/Common/config.h
    HINTS ${BBP_SDK_ROOT}/include)
endif()

set(BBP_SDK_INCLUDE_DIRS ${BBP_SDK_INCLUDE_DIR})
list(APPEND BBP_SDK_INCLUDE_DIRS "${BBP_SDK_INCLUDE_DIR}/BBP")

# Finding out the BBP-SDK version
if(BBP_SDK_INCLUDE_DIR AND EXISTS "${BBP_SDK_INCLUDE_DIR}/BBP/Common/config.h")
  set(BBP_SDK_VERSION_FILE "${BBP_SDK_INCLUDE_DIR}/BBP/Common/config.h")
endif()

if(BBP_SDK_VERSION_FILE)
  file(READ "${BBP_SDK_VERSION_FILE}" BBP_SDK_VERSION_CONTENTS)
  string(REGEX MATCH "define BBP_SDK_VERSION_MAJOR_REVISION[ \t]+([0-9]+)"
    BBP_SDK_VERSION_MAJOR ${BBP_SDK_VERSION_CONTENTS})
  string(REGEX REPLACE ".*[ \t]+([0-9]+).*" "\\1" BBP_SDK_VERSION_MAJOR
    ${BBP_SDK_VERSION_MAJOR})

  string(REGEX MATCH "define BBP_SDK_VERSION_MINOR_REVISION[ \t]+([0-9]+)"
    BBP_SDK_VERSION_MINOR ${BBP_SDK_VERSION_CONTENTS})
  string(REGEX REPLACE ".*[ \t]+([0-9]+).*" "\\1" BBP_SDK_VERSION_MINOR
    ${BBP_SDK_VERSION_MINOR})

  string(REGEX MATCH "define BBP_SDK_VERSION_PATCHLEVEL_REVISION[ \t]+([0-9]+)"
    BBP_SDK_VERSION_PATCH ${BBP_SDK_VERSION_CONTENTS})
  string(REGEX REPLACE ".*[ \t]+([0-9]+).*" "\\1" BBP_SDK_VERSION_PATCH
    ${BBP_SDK_VERSION_PATCH})

  set(BBP_SDK_VERSION "${BBP_SDK_VERSION_MAJOR}.${BBP_SDK_VERSION_MINOR}.${BBP_SDK_VERSION_PATCH}"
    CACHE INTERNAL "The version of BBP-SDK which was detected")

  if(BBP-SDK_FIND_VERSION VERSION_LESS 0.14)
    string(REGEX MATCH "define BBP_SDK_SOVERSION[ \t]+([0-9]+)"
      BBP-SDK_FIND_SOVERSION ${BBP_SDK_VERSION_CONTENTS})
    string(REGEX REPLACE ".*[ \t]+([0-9]+).*" "\\1" BBP-SDK_FIND_SOVERSION
      ${BBP-SDK_FIND_SOVERSION})
  endif()
else()
  set(BBP_SDK_FAIL TRUE)
  message(${BBP_SDK_VERSION_OUTPUT_TYPE}
    "ERROR: Can't find BBP-SDK header file config.h in ${BBP_SDK_INCLUDE_DIR}/Common/")
endif()

if(BBP-SDK_FIND_VERSION AND BBP_SDK_VERSION)
  if(BBP-SDK_FIND_VERSION_EXACT)
    if(NOT BBP_SDK_VERSION VERSION_EQUAL ${BBP-SDK_FIND_VERSION})
      set(BBP_SDK_VERSION_NOT_EXACT TRUE)
      set(BBP_SDK_FAIL TRUE)
      message(${BBP_SDK_VERSION_OUTPUT_TYPE}
        "ERROR: Version ${BBP-SDK_FIND_VERSION} of BBP-SDK is required exactly. "
        "Version ${BBP_SDK_VERSION} was found.")
    endif()
  else()
    # version is too low
    if(BBP_SDK_VERSION VERSION_LESS ${BBP-SDK_FIND_VERSION})
      set(BBP_SDK_FAIL TRUE)
      message(${BBP_SDK_VERSION_OUTPUT_TYPE}
        "ERROR: Version ${BBP-SDK_FIND_VERSION} or higher of BBP-SDK is required. "
        "Version ${BBP_SDK_VERSION} was found in ${BBP_SDK_VERSION_FILE}.")
    endif()
  endif()
endif()

# Ensure that the optional components are valid.
foreach(component ${BBP-SDK_FIND_COMPONENTS} )
  list(FIND BBP-SDK_VALID_COMPONENTS ${component} component_location)
  if(${component_location} EQUAL -1)
    message(${BBP_SDK_VERSION_OUTPUT_TYPE} "\"${component}\" is not a valid SDK component.")
  else()
    # The components don't perfectly match the libraries' names. We therefore update the components passed
    # by the user.
    if(component STREQUAL "analysis")
      list(APPEND BBP-SDK_FIND_LIBS_NAMES "filters")
      list(APPEND BBP-SDK_FIND_LIBS "BBP-SDK.filters")
    elseif(component STREQUAL "model")
      list(APPEND BBP-SDK_FIND_LIBS_NAMES "model")
      list(APPEND BBP-SDK_FIND_LIBS "BBP-SDK")
    else()
      list(APPEND BBP-SDK_FIND_LIBS_NAMES "${component}")
      list(APPEND BBP-SDK_FIND_LIBS "BBP-SDK.${component}")
    endif()
  endif()
endforeach()

# For each component, find the library
set(index 0)
foreach(COMPONENT ${BBP-SDK_FIND_LIBS})
  list(GET BBP-SDK_FIND_LIBS_NAMES ${index} COMPONENT_NAME)
  string(TOUPPER ${COMPONENT_NAME} UPPERCOMPONENT)
  unset(BBP_SDK_${UPPERCOMPONENT}_LIBRARY CACHE)
  find_library(BBP_SDK_${UPPERCOMPONENT}_LIBRARY
    NAMES ${COMPONENT} lib${COMPONENT}
    NO_DEFAULT_PATH PATHS ${BBP_SDK_ROOT}/lib ${BBP_SDK_ROOT}/lib/Release)
  # Keep track of the found components
  list(APPEND BBP_SDK_LIBRARIES "${BBP_SDK_${UPPERCOMPONENT}_LIBRARY}")
  list(APPEND BBP_SDK_COMPONENTS_FOUND "BBP_SDK_${UPPERCOMPONENT}_FOUND")
  find_package_handle_standard_args(BBP_SDK_${UPPERCOMPONENT} DEFAULT_MSG BBP_SDK_${UPPERCOMPONENT}_LIBRARY BBP_SDK_INCLUDE_DIR)
  math(EXPR index "${index}+1")
endforeach()


# find BBP-SDK data directory
unset(BBP_SDK_DATA_DIR CACHE)
find_path(BBP_SDK_DATA_DIR CMake/BBP-SDKCompilerFlags.cmake
  HINTS ${BBP_SDK_ROOT}/share/BBP-SDK${BBP-SDK_FIND_SOVERSION} ${BBP_SDK_ROOT})
# search again in latest version
if(NOT BBP_SDK_DATA_DIR)
  find_path(BBP_SDK_DATA_DIR CMake/BBP-SDKCompilerFlags.cmake
    HINTS ${BBP_SDK_ROOT}/share/BBP-SDK)
endif()
set(BBP_SDK_IDL_DIR ${BBP_SDK_DATA_DIR}/idl)

# handle the QUIETLY and REQUIRED arguments and set BBP_SDK_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(BBP-SDK DEFAULT_MSG ${BBP_SDK_COMPONENTS_FOUND} BBP_SDK_INCLUDE_DIR BBP_SDK_DATA_DIR)

# Corba needs extra include directory
if(BBP_SDK_CORBA_FOUND)
  list(APPEND BBP_SDK_INCLUDE_DIRS "${BBP_SDK_INCLUDE_DIR}/BBP/Corba")
endif()
                
if(BBP-SDK_FOUND)
  include(${BBP_SDK_DATA_DIR}/CMake/BBP-SDKCompilerFlags.cmake)
endif()

# Auxiliary dependencies
if(CMAKE_VERSION VERSION_GREATER 2.7)
  find_package(HDF5 COMPONENTS C CXX REQUIRED)
else()
  find_package(HDF5 REQUIRED)
endif()

if(HDF5_FOUND)
  list(APPEND BBP_SDK_INCLUDE_DIRS "${HDF5_INCLUDE_DIR}")
  list(APPEND BBP_SDK_LIBRARIES "${HDF5_LIBRARIES}")
else()
  set(BBP_SDK_FAIL TRUE)
endif()

if(BBP_SDK_CORBA_FOUND)
  find_package(OmniORB REQUIRED)
  list(APPEND BBP_SDK_INCLUDE_DIRS "${OMNIORB4_INCLUDE_DIR}")
  list(APPEND BBP_SDK_LIBRARIES "${OMNIORB4_LIBRARIES}")
endif()

# Wrap it up...
if(BBP_SDK_FAIL)
  # Zero out everything, we didn't meet version requirements
  set(BBP-SDK_FOUND FALSE)
  set(BBP_SDK_LIBRARIES)
  set(BBP_SDK_INCLUDE_DIRS)
else()
  if(BBP-SDK_FIND_SOVERSION VERSION_LESS 2)
    set(BBP_SDK_DEB_DEPENDENCIES "bbp-sdk${BBP-SDK_FIND_SOVERSION}")
  elseif(BBP-SDK_FIND_VERSION VERSION_LESS 0.14)
    set(BBP_SDK_DEB_DEPENDENCIES "bbpsdk${BBP-SDK_FIND_SOVERSION}")
  else()
    set(BBP_SDK_DEB_DEPENDENCIES "bbpsdk${BBP-SDK_SOVERSION}")
  endif()
endif()

set(BBP-SDK_INCLUDE_DIRS ${BBP_SDK_INCLUDE_DIRS}) # Backwards compatibility name
set(BBP-SDK_LIB ${BBP_SDK_LIBRARY}) # Backwards compatibility name
set(BBP_SDK_LIBRARY ${BBP_SDK_LIBRARIES}) # Backwards compatibility name

if(BBP-SDK_FOUND)
  message(STATUS "Found BBP-SDK ${BBP_SDK_VERSION} in ${BBP_SDK_INCLUDE_DIRS}:"
    "${BBP_SDK_LIBRARIES}")
endif()