cmake_minimum_required(VERSION 3.14)


set (PROJECT Clip)
project(${PROJECT})





include_directories(external/glm)

file(GLOB_RECURSE SRC_FILES src/*.cc)
add_executable(${PROJECT_NAME} ${SRC_FILES})

