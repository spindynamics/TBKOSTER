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

