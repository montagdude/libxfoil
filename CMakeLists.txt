# Set up project
cmake_minimum_required(VERSION 2.8.8)
project(libxfoil C)
enable_language(Fortran)
set(LIBXFOIL_VERSION 1.0)

# Determine compiler
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("Fortran compiler: " ${CMAKE_Fortran_COMPILER})

# Default compiler flags
if (Fortran_COMPILER_NAME MATCHES "gfortran")
  set (Fortran_REAL_FLAGS "-fdefault-real-8 -std=legacy")
  set (Fortran_DBLE_FLAGS "-std=legacy")
elseif (Fortran_COMPILER_NAME MATCHES "ifort")
  set (Fortran_REAL_FLAGS "-r8")
  set (Fortran_DBLE_FLAGS "")
else (Fortran_COMPILER_NAME MATCHES "gfortran")
  message(FATAL_ERROR "Fortran compiler not supported.")
endif (Fortran_COMPILER_NAME MATCHES "gfortran")

set(CMAKE_C_FLAGS_DEBUG "-g -Wall -fPIC")
set(CMAKE_C_FLAGS_RELEASE "-fPIC -O2")

# Build shared or static library
set (LIBRARY_TYPE "shared"
     CACHE STRING "Whether to build a shared or static library"
)

# Option to build examples
set (BUILD_EXAMPLES TRUE
     CACHE BOOL "Whether to build and install example programs."
)

# Set source files (exclude examples here)
file(GLOB SOURCESC "src/*.c")
file(GLOB SOURCESDBLE "src/*.f90")
file(GLOB SOURCESREAL "src/*.f")

# Source file flags
set_source_files_properties(${SOURCESREAL} PROPERTIES COMPILE_FLAGS ${Fortran_REAL_FLAGS})
set_source_files_properties(${SOURCESDBLE} PROPERTIES COMPILE_FLAGS ${Fortran_DBLE_FLAGS})

# include directory (needed by examples)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Build libraries
if (LIBRARY_TYPE MATCHES "static")
  add_library(xfoil STATIC ${SOURCESC} ${SOURCESDBLE} ${SOURCESREAL})
else (LIBRARY_TYPE MATCHES "static")
  add_library(xfoil SHARED ${SOURCESC} ${SOURCESDBLE} ${SOURCESREAL})
endif (LIBRARY_TYPE MATCHES "static")

# Build examples
if (BUILD_EXAMPLES)
  add_executable(fortran_example "examples/fortran_example.f90")
  target_link_libraries(fortran_example xfoil)
  set_target_properties(fortran_example PROPERTIES
                        INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
  
  add_executable(c_example "examples/c_example.c")
  target_link_libraries(c_example xfoil)
  set_target_properties(c_example PROPERTIES
                        INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
endif (BUILD_EXAMPLES)

# Set up distutils script for Python
configure_file(setup.py.in ${CMAKE_SOURCE_DIR}/setup.py)

# Installation
install(TARGETS xfoil DESTINATION lib${LIB_SUFFIX})
install(FILES "src/libxfoil.h" DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/libxfoil.mod DESTINATION include)
if (BUILD_EXAMPLES)
  install(TARGETS fortran_example DESTINATION ${CMAKE_SOURCE_DIR}/examples)
  install(TARGETS c_example DESTINATION ${CMAKE_SOURCE_DIR}/examples)
endif (BUILD_EXAMPLES)
