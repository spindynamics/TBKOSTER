add_executable (math_test main.f90)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if(ALTERNATE) #MATHIEU route
  message(STATUS "MATHIEU ROUTE")

  # LAPACK/BLAS and OpenMP libraries
  find_package(MKL)
  if(MKL_FOUND)
    message("-- LAPACK libraries: ${MKL_LAPACK} ${MKL_BLAS}")
    message("-- OpenMP libraries: ${INTEL_OMP}")
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
      set(OpenMP_Fortran_FLAGS -fopenmp)
    elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
      set(OpenMP_Fortran_FLAGS -qopenmp)
    endif()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DOpenMP_Fortran_FOUND")
    # See https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8 -m64")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
    include_directories(${MKL_INCLUDE_DIRS})

	#target_link_libraries(${PROJECT_NAME}.x ${MKL_LIBRARIES} ${CMAKE_DL_LIBS})
    target_link_libraries(math_test ${PROJECT_NAME} ${MKL_LIBRARIES} ${CMAKE_DL_LIBS})
  else()
    find_package(LAPACK REQUIRED)
    message("-- LAPACK libraries: ${LAPACK_LIBRARIES}")
    find_package(OpenMP REQUIRED)
    message("-- OpenMP libraries: ${OpenMP_Fortran_LIBRARIES}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")

    #target_link_libraries(${PROJECT_NAME}.x ${LAPACK_LIBRARIES} ${OpenMP_Fortran_LIBRARIES} ${CMAKE_DL_LIBS})
    target_link_libraries(math_test ${PROJECT_NAME} ${LAPACK_LIBRARIES} ${OpenMP_Fortran_LIBRARIES} ${CMAKE_DL_LIBS})
  endif()
  
else() # PASCAL route
  message(STATUS "PASCAL ROUTE")
  find_package(OpenMP)
  message(STATUS "OpenMP libraries: ${OpenMP_Fortran_LIBRARIES}")
  if (OpenMP_Fortran_FOUND)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DOpenMP_Fortran_FOUND")
  endif()

  set(BLA_F95 OFF) #put on if Intel LAPACK 95 interface is needed

  find_package(BLAS REQUIRED)
  message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")

  find_package(LAPACK REQUIRED)
  message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")

  if(BLAS95_FOUND)
    message(STATUS "BLAS Fortran 95 interface Found")
    message(STATUS "BLAS95 libraries: ${BLAS95_LIBRARIES}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DBLAS95_FOUND")
  endif()
  if(LAPACK95_FOUND)
    message(STATUS "LAPACK Fortran 95 interface Found")
    message(STATUS "LAPACK95 libraries: ${LAPACK95_LIBRARIES}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DLAPACK95_FOUND")
  endif()

  # MKL machine dependant stuffs  
  # check https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor to comply with
  if (BLAS_FOUND AND BLAS95_FOUND)
    if (BLA_VENDOR MATCHES "^Intel")
      if(NOT DEFINED ENV{MKL_IL})
        set(ENV{MKL_IL} ilp64)
      elseif(NOT $ENV{MKL_IL} STREQUAL lp64 OR NOT $ENV{MKL_IL} STREQUAL ilp64)
        message(FATAL_ERROR "Unvalid option for MKL_IL")
      endif()

      # test existance of directories
      if(NOT EXISTS "$ENV{MKLROOT}/include")
        set(ENV{MKLROOT} "$ENV{MKLROOT}/mkl")
      endif()

      # 32-bits integer
      if ($ENV{MKL_IL} STREQUAL lp64)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -I$ENV{MKLROOT}/include -I$ENV{MKLROOT}/include/intel64/lp64")
      endif()

      if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      # 64-bits integer
        if ($ENV{MKL_IL} STREQUAL ilp64)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -fdefault-integer-8 -I$ENV{MKLROOT}/include -I$ENV{MKLROOT}/include/intel64/ilp64")
        endif()
      endif()

      if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      # 64-bits integer
        if ($ENV{MKL_IL} STREQUAL ilp64)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -i8 -I$ENV{MKLROOT}/include -I$ENV{MKLROOT}/include/intel64/ilp64")
        endif()
        if (APPLE)
          set(MKL_LIBRARY_DIRS "$ENV{MKLROOT}/lib")
          find_library(MKL_CORE      mkl_core.a               HINTS ${MKL_LIBRARY_DIRS})
          find_library(MKL_COMPILER  mkl_intel_$ENV{MKL_IL}.a HINTS ${MKL_LIBRARY_DIRS})
          find_library(MKL_THREADING mkl_intel_thread.a       HINTS ${MKL_LIBRARY_DIRS})
        endif()
        if (UNIX AND NOT APPLE)
          set(MKL_LIBRARY_DIRS "$ENV{MKLROOT}/lib/intel64")
          find_library(MKL_CORE      libmkl_core.a               HINTS ${MKL_LIBRARY_DIRS})
          find_library(MKL_COMPILER  libmkl_intel_$ENV{MKL_IL}.a HINTS ${MKL_LIBRARY_DIRS})
          find_library(MKL_THREADING libmkl_intel_thread.a       HINTS ${MKL_LIBRARY_DIRS})
          set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-intel -Wl,--start-group ${MKL_CORE} ${MKL_COMPILER} ${MKL_THREADING} -Wl,--end-group")
        endif()
      endif()
    endif()
  endif()
  target_link_libraries(math_test ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} OpenMP::OpenMP_Fortran Threads::Threads ${CMAKE_DL_LIBS})
endif(ALTERNATE)