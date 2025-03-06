cmake_minimum_required(VERSION 3.14)
project(Bridger)

set(LOGGING ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



# Check for GMSH presence first
set(GMSH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gmsh)
if(NOT EXISTS ${GMSH_DIR})
  if(WIN32)
    message(FATAL_ERROR 
      "ERROR: GMSH directory not found at ${GMSH_DIR}.\n"
      "Please download GMSH from https://gmsh.info/bin/Windows/gmsh-4.13.1-Windows64-sdk.zip\n"
      "and extract it at the project root directory before running CMake again.")
  else()
    # For Mac/Linux we'll download automatically
    file(MAKE_DIRECTORY ${GMSH_DIR})
    
    # Determine platform and download appropriate SDK
    if(APPLE)
      execute_process(COMMAND uname -m OUTPUT_VARIABLE ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
      if(ARCH STREQUAL "arm64")
        set(GMSH_URL "https://gmsh.info/bin/macOS/gmsh-4.13.1-MacOSARM-sdk.tgz")
      else()
        set(GMSH_URL "https://gmsh.info/bin/macOS/gmsh-4.13.1-MacOSX-sdk.tgz")
      endif()
      set(GMSH_ARCHIVE "${GMSH_DIR}/gmsh.tgz")
    elseif(UNIX)
      set(GMSH_URL "https://gmsh.info/bin/Linux/gmsh-4.13.1-Linux64-sdk.tgz")
      set(GMSH_ARCHIVE "${GMSH_DIR}/gmsh.tgz")
    endif()
    
    # Download GMSH
    message(STATUS "Downloading GMSH from ${GMSH_URL}")
    file(DOWNLOAD ${GMSH_URL} ${GMSH_ARCHIVE} SHOW_PROGRESS)
    
    # Extract archive
    message(STATUS "Extracting GMSH...")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${GMSH_ARCHIVE} WORKING_DIRECTORY ${GMSH_DIR})
    
    # Remove downloaded archive after extraction
    file(REMOVE ${GMSH_ARCHIVE})
  endif()
endif()

# Find extracted GMSH directory
file(GLOB GMSH_SDK_DIR "${GMSH_DIR}/gmsh-*")
set(GMSH_INCLUDE_DIR "${GMSH_SDK_DIR}/include")
set(GMSH_LIB_DIR "${GMSH_SDK_DIR}/lib")



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

# Adding solver
add_subdirectory(solver)


# PARCE QUE C'EST NOTRE PROJEEET
file(GLOB_RECURSE SRC_FILES src/*.c)
add_executable(bridge ${SRC_FILES})
target_link_libraries(bridge PRIVATE glad glfw log solver)
if (LOGGING)
  target_compile_definitions(bridge PRIVATE SHOULD_LOG)
endif()

if(WIN32)
  target_link_libraries(bridge PRIVATE ${GMSH_LIB_DIR}/gmsh.lib)
else()
  target_link_libraries(bridge PRIVATE ${GMSH_LIB_DIR}/libgmsh${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()
target_include_directories(bridge PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(bridge PRIVATE solver/include)
