#/**************************************************************************
 #  CMakeLists.txt  --  This file is part of DFTCXX.                     *
 #                                                                       *
 #  Copyright (C) 2016, Ivo Filot                                        *
 #                                                                       *
 #  DFTCXX is free software:                                             *
 #  you can redistribute it and/or modify it under the terms of the      *
 #  GNU General Public License as published by the Free Software         *
 #  Foundation, either version 3 of the License, or (at your option)     *
 #  any later version.                                                   *
 #                                                                       *
 #  DFTCXX is distributed in the hope that it will be useful,            *
 #  but WITHOUT ANY WARRANTY; without even the implied warranty          *
 #  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 #  See the GNU General Public License for more details.                 *
 #                                                                       *
 #  You should have received a copy of the GNU General Public License    *
 #  along with this program.  If not, see http://www.gnu.org/licenses/.  *
 #                                                                       *
#**************************************************************************/

# set minimum cmake requirements
cmake_minimum_required(VERSION 2.8)
project (dftcxx)

# add custom directory to look for .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules )

# prepare configuration file
SET(VERSION_MAJOR "1")
SET(VERSION_MINOR "2")
SET(VERSION_MICRO "0")
configure_file(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_SOURCE_DIR}/config.h @ONLY)

# Enable release build
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Check for OpenMP
if(NOT WIN32)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        option(HAS_OPENMP "OpenMP enabled" ON)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
else()
    set(HAS_OPENMP OFF)
endif()

# set Boost
if(APPLE)
    SET(BOOST_INCLUDEDIR "/usr/local/include")
    SET(BOOST_LIBRARYDIR "/usr/local/lib")
ELSE()
    SET(BOOST_INCLUDEDIR "/usr/include")
    SET(BOOST_LIBRARYDIR "/usr/lib/x86_64-linux-gnu")
ENDIF()
SET (Boost_NO_SYSTEM_PATHS ON)
SET (Boost_USE_MULTITHREADED ON)
SET (Boost_USE_STATIC_LIBS ON)
SET (Boost_USE_STATIC_RUNTIME OFF)
SET (BOOST_ALL_DYN_LINK OFF)

find_package(Boost COMPONENTS regex iostreams filesystem REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(CPPUNIT REQUIRED) # for unit tests
pkg_check_modules(EIGEN eigen3 REQUIRED)
pkg_check_modules(TCLAP tclap REQUIRED)

# Set include folders
include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_BINARY_DIR}
                    ${TCLAP_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIRS}
                    ${EIGEN_INCLUDE_DIRS}
                    ${CPPUNIT_INCLUDE_DIR})

# use C++17 and use native compilation
add_definitions(-std=c++17)
add_definitions(-march=native)
if(OPENMP_FOUND)
    option(HAS_OPENMP "OpenMP enabled" ON)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Add sources
file(GLOB SOURCES "*.cpp")
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/dftcxx.cpp)
add_library(dftcxxsources STATIC ${SOURCES})
add_executable(dftcxx ${CMAKE_CURRENT_SOURCE_DIR}/dftcxx.cpp)

# Link libraries
target_link_libraries(dftcxx dftcxxsources ${Boost_LIBRARIES})

if(NOT DISABLE_TEST)
    message("[USER] Performing unit tests")
    enable_testing ()
    add_subdirectory("test")
else()
    message("[USER] Testing routine disabled")
endif()

###
# Installing
##
install (TARGETS dftcxx DESTINATION bin)
