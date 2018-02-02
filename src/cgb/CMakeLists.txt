include_directories(${EXTERNAL_INCLUDES})

add_library(cgb SHARED
  cgb.cc
  ${LIBTRIANGLE_OBJECT}
)

set_source_files_properties(
  ${LIBTRIANGLE_OBJECT} PROPERTIES
  EXTERNAL_OBJECT TRUE
  GENERATED TRUE
)
