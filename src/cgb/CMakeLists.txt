include_directories(${EXTERNAL_INCLUDES} ../../dependencies/eigen)
link_directories(${EXTERNAL_LIBRARIES})

add_library(cgb SHARED
  bezier.cc
  cgb.cc
  curved-domain.cc
  harmonic.cc
  lsq-plane.cc
  ${LIBTRIANGLE_OBJECT}
  )

set_source_files_properties(
  ${LIBTRIANGLE_OBJECT} PROPERTIES
  EXTERNAL_OBJECT TRUE
  GENERATED TRUE
  )

target_link_libraries(cgb
  ${LIBGEOM_LIBRARIES}
  ${LIBMEC_LIBRARIES}
  )
