cmake_minimum_required(VERSION 3.0.0)

project(fccphat)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PC_FFTW3 REQUIRED fftw3f)

include_directories("${PROJECT_SOURCE_DIR}/include")

# Compile in release mode if not set
if(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

option(USE_SIMD "Use SIMD for fcc (the compiler must support vector extensions)" OFF)
if (USE_SIMD)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    ADD_DEFINITIONS(-DFCCPHAT_USE_SIMD)
  else()
    message(FATAL_ERROR "The option USE_SIMD is not compatible with ${CMAKE_CXX_COMPILER_ID}")
  endif()
endif(USE_SIMD)

set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-march=native -Wall -ffast-math")

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")

set(SRC_FCC
  src/signal.c
  src/system.c
)

add_library(fccphat SHARED
  ${SRC_FCC}
)

target_link_libraries(fccphat
  ${PC_FFTW3_LIBRARIES}
  m
)

add_executable(tdoa
  demo/tdoa.c
)
target_link_libraries(tdoa
  fccphat
)

add_executable(speed
  demo/speed.c
)
target_link_libraries(speed
  fccphat
)

option(FCC_ENABLE_TESTS "Build tests" OFF)
if (FCC_ENABLE_TESTS)
  add_subdirectory(test)
endif (FCC_ENABLE_TESTS)

