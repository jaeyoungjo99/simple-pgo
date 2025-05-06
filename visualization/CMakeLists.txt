add_library(vis_lib
    viewer.c
    g2o_parser.c
)

target_include_directories(vis_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/core
    ${CMAKE_SOURCE_DIR}/types
    ${CMAKE_SOURCE_DIR}/visualization
)

# 이 부분 추가: OpenGL 관련 직접 링크
target_link_libraries(vis_lib
    ${GLFW_LIBRARIES}
    ${OPENGL_LIBRARIES}
    glut
    GLU
    GL
    m
)
