include_directories(${EXTERNAL_INCLUDES} ../cgb)

add_executable(example example.cc)

target_link_libraries(example
  ${LIBGEOM_LIBRARIES}
  ${LIBHARMONIC_LIBRARIES}
  cgb)
