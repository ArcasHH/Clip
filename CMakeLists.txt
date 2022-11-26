cmake_minimum_required(VERSION 3.14)


set (PROJECT Clip)
project(${PROJECT})

file(GLOB_RECURSE SRC_FILES src/*.cc)
add_executable(${PROJECT_NAME}  external/glad/src/glad.c  ${SRC_FILES} ${SRC_FILES_EXT})

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

add_subdirectory(external/glfw)
target_link_libraries(${PROJECT_NAME} glfw)

include_directories(external/glad/include)

include_directories(external/glm)
include_directories(include)



