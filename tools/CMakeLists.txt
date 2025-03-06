if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fcheck=all -fbacktrace -DDEBUG")
  set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG -Ofast")
#  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -fcheck=all -fbacktrace")
endif()
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DINTEL_COMPILER=true")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -check bounds -traceback")
endif()

set(BLA_F95 OFF) #put on if Intel LAPACK 95 interface is needed

find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
  set_config_option(HAVE_BLAS "BLAS")
endif()

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
  set_config_option(HAVE_LAPACK "LAPACK")
endif()

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

# Add bands executable
add_executable(
  bands.x
  bands.f90
)

# Add pdos executable
add_executable(
  pdos.x
  pdos.f90
)

# Add mae executable
add_executable(
  mft.x
  mft.f90
)

# Add Fermi surface executable
add_executable(
  fermisurface.x
  fermisurface.f90
)

# Add build kmesh executable
add_executable(
  build_kpoints.x
  build_kpoints.f90
)

# Add build kmesh executable
add_executable(
  dat2gnu.x
  dat2gnu.f90
)

# Add build kmesh executable
add_executable(
  col2nc.x
  col2nc.f90
)

target_link_libraries(bands.x ${PROJECT_NAME} ${CMAKE_DL_LIBS})
target_link_libraries(pdos.x ${PROJECT_NAME} ${CMAKE_DL_LIBS})
target_link_libraries(mft.x ${PROJECT_NAME} ${CMAKE_DL_LIBS})
target_link_libraries(fermisurface.x ${PROJECT_NAME} ${CMAKE_DL_LIBS})
target_link_libraries(build_kpoints.x ${PROJECT_NAME} ${CMAKE_DL_LIBS})
target_link_libraries(dat2gnu.x ${PROJECT_NAME} ${CMAKE_DL_LIBS})
target_link_libraries(col2nc.x ${PROJECT_NAME} ${BLAS95_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES} ${CMAKE_DL_LIBS})
