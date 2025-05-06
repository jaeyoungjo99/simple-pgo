add_library(core_lib
    graph.c
    optimizer.c
)

target_include_directories(core_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)
