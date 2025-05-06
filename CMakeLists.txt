cmake_minimum_required(VERSION 3.10)
project(simple_pgo C)

set(CMAKE_C_STANDARD 99)

# GLFW & OpenGL for visualization
find_package(OpenGL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(GLFW REQUIRED glfw3)

include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${GLFW_INCLUDE_DIRS}
)

add_subdirectory(core)
add_subdirectory(types)
add_subdirectory(visualization)

add_executable(simple_pgo main.c)
target_link_libraries(simple_pgo
    core_lib
    types_lib
    vis_lib
    ${GLFW_LIBRARIES}
    ${OPENGL_LIBRARIES}
)
