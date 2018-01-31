# Set up project
cmake_minimum_required(VERSION 2.8.8)
project(libxfoil)
enable_language(Fortran)
set(LIBXFOIL_VERSION 0.1)

# Determine compiler
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("Fortran compiler: " ${CMAKE_Fortran_COMPILER})

# Enable OpenMP by default for Release builds, but disable by default for Debug
if (CMAKE_BUILD_TYPE MATCHES "Release")
  set (OPENMP_FLAG "-fopenmp")
  if (NOT DEFINED ENABLE_OPENMP)
    set(ENABLE_OPENMP TRUE CACHE BOOL "Whether to build with OpenMP support.")
  endif (NOT DEFINED ENABLE_OPENMP)
  if (NOT ENABLE_OPENMP)
    MESSAGE(WARNING "Disabling OpenMP support since ENABLE_OPENMP=FALSE")
    set (OPENMP_FLAG "")
  endif (NOT ENABLE_OPENMP)
elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
  set (OPENMP_FLAG "")
  if (NOT DEFINED ENABLE_OPENMP)
    set(ENABLE_OPENMP FALSE CACHE BOOL "Whether to build with OpenMP support.")
  endif (NOT DEFINED ENABLE_OPENMP)
  if (ENABLE_OPENMP)
    MESSAGE(WARNING "Enabling OpenMP support for debug release since ENABLE_OPENMP=TRUE")
    set (OPENMP_FLAG "-fopenmp")
  endif (ENABLE_OPENMP)
endif (CMAKE_BUILD_TYPE MATCHES "Release")

# Default compiler flags
if (Fortran_COMPILER_NAME MATCHES "gfortran")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fPIC ${OPENMP_FLAG}")
  set (CMAKE_Fortran_FLAGS_DEBUG "-fPIC -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall ${OPENMP_FLAG}")
  set (Fortran_REAL8_FLAG "-fdefault-real-8")
elseif (Fortran_COMPILER_NAME MATCHES "ifort")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${OPENMP_FLAG}")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -warn all ${OPENMP_FLAG}")
  set (Fortran_REAL8_FLAG "-r8")
else (Fortran_COMPILER_NAME MATCHES "gfortran")
  message(FATAL_ERROR "Fortran compiler not supported.")
endif (Fortran_COMPILER_NAME MATCHES "gfortran")

set(CMAKE_C_FLAGS_DEBUG "-g -Wall -fPIC")
set(CMAKE_C_FLAGS_RELEASE "-fPIC -O2")

# Option to build examples
set (BUILD_EXAMPLES TRUE
     CACHE BOOL "Whether to build and install example programs."
)

# Set source files (exclude examples here)
file(GLOB SOURCESDBLE "src/*.f90")
file(GLOB SOURCESREAL "src/*.f")
list(REMOVE_ITEM SOURCESDBLE "${CMAKE_SOURCE_DIR}/src/xfoil_interface_wrap.f90")

# Real-8 flag for xfoil and naca456 source files
set_source_files_properties(${SOURCESREAL} PROPERTIES COMPILE_FLAGS ${Fortran_REAL8_FLAG})

# include directory (needed by examples)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Build shared library
add_library("xfoil" SHARED ${SOURCESREAL} ${SOURCESDBLE})

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
install(FILES "src/xfoil_interface.h" "src/xfoil_interface_wrap.f90" DESTINATION include)
if (BUILD_EXAMPLES)
  install(TARGETS fortran_example DESTINATION ${CMAKE_SOURCE_DIR}/examples)
  install(TARGETS c_example DESTINATION ${CMAKE_SOURCE_DIR}/examples)
endif (BUILD_EXAMPLES)
