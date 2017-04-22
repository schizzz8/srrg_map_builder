add_executable(merging_test merging_test.cpp)
target_link_libraries(merging_test 
  srrg_map_builder_library
  ${catkin_LIBRARIES}
)

add_executable(refining_test refining_test.cpp)
target_link_libraries(refining_test 
  srrg_map_builder_library
  ${catkin_LIBRARIES}
)
