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
  set (Fortran_REAL8_FLAG "-fdefault-real-8")
elseif (Fortran_COMPILER_NAME MATCHES "ifort")
  set (Fortran_REAL8_FLAG "-r8")
else (Fortran_COMPILER_NAME MATCHES "gfortran")
  message(FATAL_ERROR "Fortran compiler not supported.")
endif (Fortran_COMPILER_NAME MATCHES "gfortran")

set(CMAKE_C_FLAGS_DEBUG "-g -Wall -fPIC")
set(CMAKE_C_FLAGS_RELEASE "-fPIC -O2")

# Build shared or static library
set (SHARED_LIBRARY TRUE
     CACHE BOOL "Whether to build a shared library")
set (STATIC_LIBRARY FALSE
     CACHE BOOL "Whether to build a static library")


# Option to build examples
set (BUILD_EXAMPLES TRUE
     CACHE BOOL "Whether to build and install example programs."
)

# Set source files (exclude examples here)
file(GLOB SOURCESC "src/*.c")
file(GLOB SOURCESDBLE "src/*.f90")
file(GLOB SOURCESREAL "src/*.f")

# Real-8 flag for xfoil source files
set_source_files_properties(${SOURCESREAL} PROPERTIES COMPILE_FLAGS ${Fortran_REAL8_FLAG})

# include directory (needed by examples)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Build libraries
if (SHARED_LIBRARY)
  add_library("xfoil_shared" SHARED ${SOURCESC} ${SOURCESDBLE} ${SOURCESREAL})
  set_target_properties(xfoil_shared PROPERTIES OUTPUT_NAME xfoil)
endif (SHARED_LIBRARY)
if (STATIC_LIBRARY)
  add_library("xfoil_static" STATIC ${SOURCESC} ${SOURCESDBLE} ${SOURCESREAL})
  set_target_properties(xfoil_static PROPERTIES OUTPUT_NAME xfoil)
endif (STATIC_LIBRARY)

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
if (SHARED_LIBRARY)
  install(TARGETS xfoil_shared DESTINATION lib${LIB_SUFFIX})
endif (SHARED_LIBRARY)
if (STATIC_LIBRARY)
  install(TARGETS xfoil_static DESTINATION lib${LIB_SUFFIX})
endif (STATIC_LIBRARY)
install(FILES "src/libxfoil.h" DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/libxfoil.mod DESTINATION include)
if (BUILD_EXAMPLES)
  install(TARGETS fortran_example DESTINATION ${CMAKE_SOURCE_DIR}/examples)
  install(TARGETS c_example DESTINATION ${CMAKE_SOURCE_DIR}/examples)
endif (BUILD_EXAMPLES)
