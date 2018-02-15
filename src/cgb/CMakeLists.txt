include_directories(${EXTERNAL_INCLUDES} ../../dependencies/eigen)
link_directories(${EXTERNAL_LIBRARIES})

add_library(cgb SHARED
  cgb.cc
  ${LIBTRIANGLE_OBJECT}
  )

set_source_files_properties(
  ${LIBTRIANGLE_OBJECT} PROPERTIES
  EXTERNAL_OBJECT TRUE
  GENERATED TRUE
  )

target_link_libraries(cgb
  ${LIBGEOM_LIBRARIES}
  ${LIBHARMONIC_LIBRARIES}
  )
