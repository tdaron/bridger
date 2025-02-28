cmake_minimum_required(VERSION 3.14)
project(Bridger)

set(LOGGING ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Adding GLFW
add_subdirectory(glfw)  
target_include_directories(glfw PRIVATE glfw/include)

# Adding GLAD
add_library(glad STATIC glad/src/glad.c)
target_include_directories(glad PUBLIC glad/include)

# Adding log
add_library(log STATIC log/src/log.c)
target_include_directories(log PUBLIC log/src/)
target_compile_definitions(log PRIVATE LOG_USE_COLOR)


# PARCE QUE C'EST NOTRE PROJEEET
file(GLOB_RECURSE SRC_FILES src/*.c)
add_executable(bridge ${SRC_FILES})
target_link_libraries(bridge PRIVATE glad glfw log)
if (LOGGING)
  target_compile_definitions(bridge PRIVATE SHOULD_LOG)
endif()
