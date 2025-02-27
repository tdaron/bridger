cmake_minimum_required(VERSION 3.14)
project(Bridger)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# glad target
add_library(glad STATIC glad/src/glad.c)
target_include_directories(glad PUBLIC glad/include)

# Add GLFW as a dependency
add_subdirectory(glfw)  # Replace with the correct path to GLFW source

# macOS-specific settings
if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64")  # Supports both Intel and Apple Silicon Macs
    set(CMAKE_CXX_STANDARD 17)  # C++17 is recommended for macOS
    target_link_libraries(bridge PRIVATE "-framework Cocoa" "-framework OpenGL" "-framework IOKit" "-framework CoreVideo")
endif()

# bridge target
add_executable(bridge src/main.c src/gl_init.c src/read_mesh.c)
target_link_libraries(bridge PRIVATE glad glfw)
target_include_directories(bridge PRIVATE glad/include)
target_include_directories(bridge PRIVATE glfw/include)
