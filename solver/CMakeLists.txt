cmake_minimum_required(VERSION 3.14)
project(Solver)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB_RECURSE SRC_FILES src/*.c)
list(FILTER SRC_FILES EXCLUDE REGEX ".*/main\.c$")

add_library(solver ${SRC_FILES})

file(GLOB_RECURSE ALL_SRC_FILES src/*.c)
add_executable(solver-bin ${ALL_SRC_FILES})
target_include_directories(solver-bin PRIVATE include)
