add_library(srrg_map_builder_library SHARED
	map_builder.h map_builder.cpp
)

target_link_libraries(srrg_map_builder_library
	${catkin_LIBRARIES}
)
