include_directories(${EXTERNAL_INCLUDES} ../cgb)
link_directories(${EXTERNAL_LIBRARIES})

add_executable(example example.cc)

target_link_libraries(example cgb)
