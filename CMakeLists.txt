#[[Copyright (c) 2020 Jonas Ellert

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.]]

cmake_minimum_required(VERSION 3.12.4 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(par-max-dicut)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

# Compiler flags
set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -fopenmp -Wall -Wextra -Wpedantic -fdiagnostics-color=auto")
set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb -DDEBUG")

# Libraries
set(PAR_MAX_CUT_LIBRARIES stdc++fs)

# Add ips4o
include_directories(${PROJECT_SOURCE_DIR}/external/ips4o)

# Add KaMinPar
add_subdirectory(${PROJECT_SOURCE_DIR}/external/KaMinPar)
set(kaminpar_libs kaminpar tbb pthread atomic omp numa)
set(PAR_MAX_CUT_LIBRARIES ${PAR_MAX_CUT_LIBRARIES} ${kaminpar_libs})
include_directories("${PROJECT_SOURCE_DIR}/external/KaMinPar/library")

# Add tlx command line parser
add_library(tlx_command_line
        external/tlx/tlx/cmdline_parser.cpp
        external/tlx/tlx/logger/core.cpp
        external/tlx/tlx/string/parse_si_iec_units.cpp)
set(PAR_MAX_CUT_LIBRARIES ${PAR_MAX_CUT_LIBRARIES} tlx_command_line)
include_directories("${PROJECT_SOURCE_DIR}/external/tlx")

# Add Gurobi Solver
find_package(GUROBI)
if(GUROBI_FOUND)
    include_directories(${GUROBI_INCLUDE_DIRS})
    set(PAR_MAX_CUT_LIBRARIES ${PAR_MAX_CUT_LIBRARIES} ${GUROBI_CXX_LIBRARY} ${GUROBI_CXX_DEBUG_LIBRARY} ${GUROBI_LIBRARY})
    add_definitions(-DUSE_GUROBI=${GUROBI_FOUND})
endif()

# Add Mosek Solver
find_package(MOSEK)
if(MOSEK_FOUND)
    include_directories(${MOSEK_INCLUDE_DIRS})
    message(STATUS "MOSEK_INCLUDE_DIRS = ${MOSEK_INCLUDE_DIRS}")
    message(STATUS "MOSEK_LIBRARY = ${MOSEK_LIBRARY}")
    message(STATUS "FUSION_LIBRARY = ${FUSION_LIBRARY}")
    set(PAR_MAX_CUT_LIBRARIES ${PAR_MAX_CUT_LIBRARIES} ${MOSEK_LIBRARY} ${FUSION_LIBRARY})
    add_definitions(-DUSE_MOSEK=${MOSEK_FOUND})
endif()

include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(benchmark "src/benchmark.cpp")
target_link_libraries(benchmark ${PAR_MAX_CUT_LIBRARIES})

# add library

message(STATUS "Adding library 'par-max-dicut_lib'")

add_library(par-max-dicut_lib INTERFACE)
target_include_directories(par-max-dicut_lib INTERFACE "include/")
if (GUROBI_FOUND)
    target_include_directories(par-max-dicut_lib INTERFACE ${GUROBI_INCLUDE_DIRS})
endif()
if (MOSEK_FOUND)
    target_include_directories(par-max-dicut_lib INTERFACE ${MOSEK_INCLUDE_DIRS})
endif()
target_link_libraries(par-max-dicut_lib INTERFACE ${PAR_MAX_CUT_LIBRARIES})

message(STATUS "Library 'par-max-dicut_lib' successfully added")

MESSAGE(STATUS "Built Type: " ${CMAKE_BUILD_TYPE} )

################################################################################
