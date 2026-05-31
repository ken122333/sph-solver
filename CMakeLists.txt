cmake_minimum_required(VERSION 3.20)
project(sph_solver LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# GLFW
FetchContent_Declare(glfw
  GIT_REPOSITORY https://github.com/glfw/glfw.git
  GIT_TAG 3.3.9
)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)

# Python (use launcher)
set(Python3_EXECUTABLE "py")
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# GLAD (generated)
FetchContent_Declare(glad
  GIT_REPOSITORY https://github.com/Dav1dde/glad.git
  GIT_TAG v2.0.6
)
FetchContent_MakeAvailable(glad)
add_subdirectory("${glad_SOURCE_DIR}/cmake" glad_cmake)

glad_add_library(glad
  REPRODUCIBLE
  EXCLUDE_FROM_ALL
  LOADER
  API gl:core=3.3
)

# OpenGL
find_package(OpenGL REQUIRED)

# Your executable
add_executable(sph
  src/main.cpp
  src/sph/solver.cpp
 "src/sph/particle.cpp" 
 "src/sph/spatial_hash_grid.cpp"
 "src/sph/render.cpp")

target_include_directories(sph PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${glfw_SOURCE_DIR}/include
)

target_link_libraries(sph PRIVATE
  glfw
  glad
  OpenGL::GL
)

# Inherit include directories from the glad target (this is the key)
get_target_property(_glad_includes glad INCLUDE_DIRECTORIES)
if(_glad_includes)
  target_include_directories(sph PRIVATE ${_glad_includes})
endif()
