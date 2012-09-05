#Ecole Polytechnique Federale de Lausanne
#Brain Mind Institute,
#Blue Brain Project
#(c) 2006-2008. All rights reserved.

# _____________________________________________________________________________
#
# BBP-SDK 
# _____________________________________________________________________________
#

cmake_minimum_required(VERSION 2.8)

MESSAGE("Entering FindBBP")
set(BBP_VERSION "0.14.0")

# PATH ________________________________________________________________________

find_path(BBP-SDK_PATH include/BBP/common.h
    ${CMAKE_SOURCE_DIR}/../../../BBP-SDK/
    ${CMAKE_SOURCE_DIR}/../../BBP-SDK/
    "C:/Program Files/BBP-SDK"
    /opt
    /sw/BBP-SDK
    /opt/BBP-SDK
    D:/cmakebuild/buildyard/install
)

if (BBP-SDK_PATH)
    set (BBP-SDK_FOUND TRUE)
endif (BBP-SDK_PATH)

# HEADERS _____________________________________________________________________

if (BBP-SDK_FOUND)
    set (BBP-SDK_INCLUDE_DIRS ${BBP-SDK_PATH}/include)
    mark_as_advanced (BBP-SDK_INCLUDE_DIR)


# DYNAMIC OR STATIC LIBRARY ___________________________________________________

    find_library(BBP-SDK_LIB 
        NAMES BBP-SDK.${BBP_VERSION} 
        PATHS ${BBP-SDK_PATH}/lib
              ${BBP-SDK_PATH}/lib/Release
    )
    find_library(BBP-SDK_LIB_DEBUG
        NAMES BBP-SDK.dbg
        PATHS ${BBP-SDK_PATH}/lib
              ${BBP-SDK_PATH}/lib/Debug
    )
    mark_as_advanced(BBP-SDK_LIB)  
    mark_as_advanced(BBP-SDK_LIB_DEBUG)    

    find_library(BBP-SDK_CORBA_LIB 
        NAMES BBP-SDK-CORBA
        PATHS ${BBP-SDK_PATH}/lib
              ${BBP-SDK_PATH}/lib/Release
    )
    find_library(BBP-SDK_CORBA_LIB_DEBUG
        NAMES BBP-SDK-CORBA.dbg
        PATHS ${BBP-SDK_PATH}/lib
              ${BBP-SDK_PATH}/lib/Debug
    )
    mark_as_advanced(BBP-SDK_CORBA_LIB)  
    mark_as_advanced(BBP-SDK_CORBA_LIB_DEBUG)

# COMPILE OPTIONS _____________________________________________________________
    if (EXISTS ${BBP-SDK_PATH}/lib/BBP-SDK.Debug.cmake AND
        "${CMAKE_BUILD_TYPE}" MATCHES "Debug")
        include(${BBP-SDK_PATH}/lib/BBP-SDK.Debug.cmake)
        add_definitions(${BBP-SDK_DEFINITIONS})
    endif (EXISTS ${BBP-SDK_PATH}/lib/BBP-SDK.Debug.cmake AND
           "${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    if (EXISTS ${BBP-SDK_PATH}/lib/BBP-SDK.Release.cmake AND
        "${CMAKE_BUILD_TYPE}" MATCHES "Release|^$")
        include(${BBP-SDK_PATH}/lib/BBP-SDK.Release.cmake)
        add_definitions(${BBP-SDK_DEFINITIONS})
    endif (EXISTS ${BBP-SDK_PATH}/lib/BBP-SDK.Release.cmake AND
         "${CMAKE_BUILD_TYPE}" MATCHES "Release|^$")
    if (EXISTS ${BBP-SDK_PATH}/lib/BBP-SDK.RelWithDebInfo.cmake AND
        "${CMAKE_BUILD_TYPE}" MATCHES "RelWithDebInfo")
        include(${BBP-SDK_PATH}/lib/BBP-SDK.RelWithDebInfo.cmake)
        add_definitions(${BBP-SDK_DEFINITIONS})
    endif (EXISTS ${BBP-SDK_PATH}/lib/BBP-SDK.RelWithDebInfo.cmake AND
           "${CMAKE_BUILD_TYPE}" MATCHES "RelWithDebInfo")

    set(BBP-SDK_INCLUDE_DIRS ${BBP-SDK_INCLUDE_DIRS}
                             ${HDF5_INCLUDE_DIR}
                             ${XML2_INCLUDE_DIR}
                             ${Boost_INCLUDE_DIR})

# FOUND _______________________________________________________________________

   if (NOT BBP-SDK_FIND_QUIETLY)
      message(STATUS "Found BBP-SDK: ${BBP-SDK_PATH}")
   endif (NOT BBP-SDK_FIND_QUIETLY)
else (BBP-SDK_FOUND)
   if (BBP-SDK_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find BBP-SDK")
   endif (BBP-SDK_FIND_REQUIRED)
endif (BBP-SDK_FOUND)
