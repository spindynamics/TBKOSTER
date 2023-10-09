if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fcheck=all -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG -Ofast")
#  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -fcheck=all -fbacktrace")
endif()
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DINTEL_COMPILER=true")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -check bounds -traceback")
endif()

# Add bands executable
add_executable(
  bands.x
  bands.f90
)

# Add bands executable
add_executable(
  pdos.x
  pdos.f90
)

target_link_libraries(bands.x ${PROJECT_NAME} ${CMAKE_DL_LIBS})
target_link_libraries(pdos.x ${PROJECT_NAME} ${CMAKE_DL_LIBS})

