cmake_minimum_required(VERSION 2.8)

project(viewer)

find_package(OpenGL REQUIRED)

if (APPLE)
  set(CMAKE_MACOSX_RPATH 1)
endif()

if (UNIX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=gnu++11")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-variable")
  # recommended but not set by default
  # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
elseif(MSVC)
  # recommended but not set by default
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -WX")
endif()

set (VIEWER_LIBS ${OPENGL_gl_LIBRARY})

# GLFW
set(GLFW_INSTALL OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL " " FORCE)
add_subdirectory(3rd_party/glfw)
list(APPEND VIEWER_LIBS glfw)
list(APPEND VIEWER_INCLUDES 3rd_party/glfw/include)

# GLEW - not needed on OS X
# we add glew source/header directly to the build, no glew library build.
if (NOT APPLE)
  add_definitions(-DGLEW_STATIC)
  list(APPEND VIEWER_INCLUDES 3rd_party/glew/include)
  list(APPEND VIEWER_SRC 3rd_party/glew/src/glew.c)
  SOURCE_GROUP(GLEW FILES 3rd_party/glew/src/glew.c)
endif()


# vecmath include directory
include_directories(vecmath/include)
add_subdirectory(vecmath)
list (APPEND VIEWER_LIBS vecmath)
list (APPEND VIEWER_INCLUDES vecmath/include)
list (APPEND VIEWER_SRC
  src/main.cpp
  src/starter3_util.cpp
  src/camera.cpp
  src/recorder.cpp
  src/vertexrecorder.cpp
  src/gl_util.cpp
)
list (APPEND VIEWER_HEADER
  src/gl.h
  src/starter3_util.h
  src/camera.h
  src/recorder.h
  src/vertexrecorder.h
  src/gl_util.h
)

add_executable(viewer ${VIEWER_SRC} ${VIEWER_HEADER})
target_include_directories(viewer PUBLIC ${VIEWER_INCLUDES})
target_link_libraries(viewer ${VIEWER_LIBS})
