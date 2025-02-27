cmake_minimum_required(VERSION 3.14)
project(Bridger)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# glad target
add_library(glad STATIC glad/src/glad.c)
target_include_directories(glad PUBLIC glad/include)

# Add GLFW as a dependency
add_subdirectory(glfw)  # Replace with the correct path to GLFW source

# bridge target
add_executable(bridge src/main.c src/gl_init.c src/read_mesh.c)
target_link_libraries(bridge PRIVATE glad glfw)
target_include_directories(bridge PRIVATE glad/include)
target_include_directories(bridge PRIVATE glfw/include)
