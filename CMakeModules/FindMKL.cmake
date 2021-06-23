# Find the Math Kernel Library (MKL) from Intel
#
# Environment variables used by this file
# INTELROOT: Intel root directory (default: /opt/intel)
# MKLROOT  : MKL root directory (default: $ENV{INTELROOT}/mkl)
# F95ROOT  : F95 root directory (default: $ENV{INTELROOT}/mkl_gf_<VERSION> with gfortran-<VERSION>
#                                         $ENV{INTELROOT}/mkl              with ifort)
#            The BLAS/LAPACK Fortran95 interface is compiler specific.
#            With gfortran, it must be compiled in a separate directory.
#            With ifort, it is provided by default in the MKL distribution.
# MKLIL    : MKL interface layer, options: lp64 (default), ilp64
#
# Copyright (C) 2017 Mathieu César <mailto:mathieu.cesar@cea.fr>

if(UNIX AND NOT APPLE)
	# Determine Fortran compiler
	if(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU AND NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
		message(FATAL_ERROR "FindMKL.cmake: Fortran compiler not supported")
	endif()

	# Initialize environment variables
	if(NOT DEFINED ENV{INTELROOT})
		set(ENV{INTELROOT} /opt/intel)
	endif()
	if(NOT DEFINED ENV{MKLROOT})
		set(ENV{MKLROOT} $ENV{INTELROOT}/mkl)
	endif()
	if(NOT DEFINED ENV{F95ROOT})
		if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
			if(EXISTS $ENV{MKLROOT}_gf_${CMAKE_Fortran_COMPILER_VERSION})
				set(ENV{F95ROOT} $ENV{MKLROOT}_gf_${CMAKE_Fortran_COMPILER_VERSION})
				add_definitions(-DBLAS95_FOUND -DLAPACK95_FOUND)
			else()
				set(ENV{F95ROOT} $ENV{MKLROOT})
			endif()
		elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
			set(ENV{F95ROOT} $ENV{MKLROOT})
			add_definitions(-DBLAS95_FOUND -DLAPACK95_FOUND)
		endif()
	endif()
	if(NOT DEFINED ENV{MKLIL})
		set(ENV{MKLIL} lp64)
	elseif(NOT $ENV{MKLIL} STREQUAL lp64 OR NOT $ENV{MKLIL} STREQUAL ilp64)
		message(FATAL_ERROR "Option for MKLIL not valid, must be one of (lp64,ilp64)")
	endif()
	message("-- ENV{INTELROOT}=$ENV{INTELROOT}")
	message("-- ENV{MKLROOT}=$ENV{MKLROOT}")
	message("-- ENV{F95ROOT}=$ENV{F95ROOT}")
	message("-- ENV{MKLIL}=$ENV{MKLIL}")

	# Set MKL_INCLUDE_DIRS
	# Undefine include directories from the cache
	unset(MKL_INCLUDE_DIR1 CACHE)
	unset(MKL_INCLUDE_DIR2 CACHE)
	unset(MKL_INCLUDE_DIRS CACHE)
	# Find include directories
	find_path(MKL_INCLUDE_DIR1 lapack95.mod PATHS $ENV{F95ROOT}/include/intel64/$ENV{MKLIL})
	find_path(MKL_INCLUDE_DIR2 lapack95.mod PATHS $ENV{MKLROOT}/include/intel64/$ENV{MKLIL})
	# Set MKL_INCLUDE_DIRS
	set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR1} ${MKL_INCLUDE_DIR2})
	message("-- MKL_INCLUDE_DIRS=${MKL_INCLUDE_DIRS}")

	# See https://software.intel.com/en-us/mkl-linux-developer-guide-listing-libraries-on-a-link-line
	# Set MKL_LIBRARIES
	# Undefine libraries from the cache
	unset(MKL_LIBRARY_DIRS CACHE)
	unset(F95_LIBRARY_DIRS CACHE)
	unset(INTEL_OMP_LIBRARY_DIRS CACHE)
	unset(MKL_BLAS CACHE)
	unset(MKL_LAPACK CACHE)
	unset(MKL_CORE CACHE)
	unset(MKL_COMPILER CACHE)
	unset(MKL_THREADING CACHE)
	unset(INTEL_OMP CACHE)
	unset(MKL_LIBRARIES CACHE)
	# Find libraries
	set(MKL_LIBRARY_DIRS $ENV{MKLROOT}/lib/intel64)
	set(F95_LIBRARY_DIRS $ENV{F95ROOT}/lib/intel64)
	set(INTEL_OMP_LIBRARY_DIRS $ENV{INTELROOT}/lib/intel64)
	if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
		find_library(MKL_BLAS     libmkl_blas95_$ENV{MKLIL}.a   PATHS ${F95_LIBRARY_DIRS})
		find_library(MKL_LAPACK   libmkl_lapack95_$ENV{MKLIL}.a PATHS ${F95_LIBRARY_DIRS})
		find_library(MKL_COMPILER libmkl_gf_$ENV{MKLIL}.a       PATHS ${MKL_LIBRARY_DIRS})
	elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
		find_library(MKL_BLAS     libmkl_blas95_$ENV{MKLIL}.a   PATHS ${MKL_LIBRARY_DIRS})
		find_library(MKL_LAPACK   libmkl_lapack95_$ENV{MKLIL}.a PATHS ${MKL_LIBRARY_DIRS})
		find_library(MKL_COMPILER libmkl_intel_$ENV{MKLIL}.a    PATHS ${MKL_LIBRARY_DIRS})
	endif()
	find_library(MKL_THREADING libmkl_intel_thread.a PATHS ${MKL_LIBRARY_DIRS})
	find_library(MKL_CORE      libmkl_core.a         PATHS ${MKL_LIBRARY_DIRS})
	find_library(INTEL_OMP     libiomp5.a            PATHS ${INTEL_OMP_LIBRARY_DIRS})
	# Set MKL_LIBRARIES
	set(MKL_LIBRARIES ${MKL_BLAS} ${MKL_LAPACK} -Wl,--start-group ${MKL_COMPILER} ${MKL_THREADING} ${MKL_CORE} -Wl,--end-group ${INTEL_OMP})
	message("-- MKL_LIBRARIES=${MKL_LIBRARIES}")

	# Validate include directories and libraries, then set MKL_FOUND
	include(FindPackageHandleStandardArgs)
	FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_INCLUDE_DIRS MKL_BLAS MKL_LAPACK MKL_CORE MKL_COMPILER MKL_THREADING)
endif()

if(APPLE)
	if(BLA_VENDOR MATCHES "^Intel")
	SET(MKL_INCLUDE_DIRS "$ENV{MKLROOT}/include")
	SET(MKL_LIBRARY_DIRS "$ENV{MKLROOT}/lib")
	SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortan_FLAGS} -I${MKL_INCLUDE_DIRS}")
	if(NOT DEFINED ENV{MKLIL})
		set(ENV(MKLIL) lp64) # 32-bits integer
	elseif(NOT $ENV{MKLIL} STREQUAL lp64 OR NOT $ENV{MKLIL} STREQUAL ilp64)
		message(FATAL_ERROR "Option for MKLIL not valid, must be one of (lp64,ilp64)")
	endif()
	unset(MKL_CORE CACHE)
	unset(MKL_COMPILER CACHE)
	unset(MKL_THREADING CACHE)
	unset(MKL_LIBOMP CACHE)
	find_library(MKL_CORE      mkl_core.a              HINTS ${MKL_LIBRARY_DIRS})
	find_library(MKL_COMPILER  mkl_intel_$ENV{MKLIL}.a HINTS ${MKL_LIBRARY_DIRS})
	find_library(MKL_THREADING mkl_intel_thread.a      HINTS ${MKL_LIBRARY_DIRS})
	find_library(MKL_LIBOMP    iomp5.a                 HINTS $ENV{INTELROOT}/compilers_and_libraries_2017.5.220/mac/compiler/lib)
	# Validate include directories and libraries, then set MKL_FOUND
	include(FindPackageHandleStandardArgs)
	FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_INCLUDE_DIRS MKL_CORE MKL_COMPILER MKL_THREADING MKL_LIBOMP)
	endif()
endif()
