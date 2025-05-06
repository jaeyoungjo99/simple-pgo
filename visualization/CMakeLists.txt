add_library(vis_lib
    viewer.c
)

target_include_directories(vis_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${GLFW_INCLUDE_DIRS}
)

target_link_libraries(vis_lib
    ${GLFW_LIBRARIES}
    ${OPENGL_LIBRARIES}
)
