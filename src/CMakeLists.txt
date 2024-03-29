cmake_minimum_required(VERSION 2.8)
project(ConcaveGB)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(LibGeom REQUIRED)
find_package(LibMec REQUIRED)
find_package(LibTriangle REQUIRED)

set(EXTERNAL_INCLUDES
  ${LIBGEOM_INCLUDE_DIRS}
  ${LIBMEC_INCLUDE_DIRS}
  ${LIBTRIANGLE_INCLUDE_DIRS}
  )

set(EXTERNAL_LIBRARIES
  ${LIBGEOM_LIB_DIRS}
  ${LIBMEC_LIB_DIRS}
  )

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(CMAKE_COMPILER_IS_GNUCXX)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-Weverything -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

add_subdirectory(cgb)
add_subdirectory(example)
