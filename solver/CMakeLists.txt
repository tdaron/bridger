cmake_minimum_required(VERSION 3.14)
project(Solver)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set degbug mode
# set(CMAKE_BUILD_TYPE Debug)
# Set release mode
set(CMAKE_BUILD_TYPE Release)

file(GLOB_RECURSE SRC_FILES src/*.c)
list(FILTER SRC_FILES EXCLUDE REGEX ".*/main\.c$")

add_library(solver ${SRC_FILES})
target_include_directories(solver PRIVATE include)

file(GLOB_RECURSE ALL_SRC_FILES src/*.c)
add_executable(solver-bin ${ALL_SRC_FILES})
target_link_libraries(solver-bin m)
target_include_directories(solver-bin PRIVATE include)
