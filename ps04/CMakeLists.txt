# TMA4280 Supercomputing, Introduction
# Martin Ukrop (martiu@stud.ntnu.no)
# problem set 04
# 2014-02-11

project(ps04)
enable_language(CXX)
cmake_minimum_required(VERSION 2.8)

option(ENABLE_MPI "Enable MPI support?" ON)
option(ENABLE_OPENMP "Enable OpenMP support?" ON)

if(ENABLE_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    list(APPEND DEFINITIONS -DHAVE_MPI)
    list(APPEND CXX_FLAGS ${MPI_CXX_FLAGS})
    list(APPEND LIBRARIES ${MPI_CXX_LIBRARIES})
    list(APPEND INCLUDE_PATHS ${MPI_CXX_INCLUDE_PATH})
  endif()
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    list(APPEND DEFINITIONS -DHAVE_OPENMP)
    list(APPEND CXX_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endif()

list(APPEND CMAKE_CXX_FLAGS ${CXX_FLAGS})
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O2")
add_definitions(${DEFINITIONS})
include_directories(${INCLUDE_PATHS})

set(SRC ps04.cpp)
set(EXECUTABLE ps04)

add_executable(${EXECUTABLE} ${SRC})
target_link_libraries(${EXECUTABLE} ${LIBRARIES})

message(STATUS "include paths=${INCLUDE_PATHS}")
message(STATUS "compile flags=${CMAKE_CXX_FLAGS}")
message(STATUS "definitions=${DEFINITIONS}")
message(STATUS "libraries=${LIBRARIES}")

set(MPIEXEC_ARGS -mca btl ^openib)

enable_testing()
# serial code
add_test(serial ${PROJECT_BINARY_DIR}/${EXECUTABLE} 0.1)
# MPI code with 1 processor
add_test(mpi-1 ${MPIEXEC} ${MPIEXEC_ARGS} ${PROJECT_BINARY_DIR}/${EXECUTABLE} 0.1)
# MPI code with 2 processors
add_test(mpi-2 ${MPIEXEC} ${MPIEXEC_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2 ${PROJECT_BINARY_DIR}/${EXECUTABLE} 0.1)

# convenience target for verbose testing
add_custom_target(verbose-test ctest -V WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
