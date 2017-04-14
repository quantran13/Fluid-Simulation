cmake_minimum_required(VERSION 3.5)
project(fluid_simulation)

set(CMAKE_C_STANDARD 11)

set(SOURCE_FILES src/main.c src/graphic.c src/fluid.c includes/fluid.h includes/graphic.h src/utility.cpp includes/utility.h)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories(includes/)
include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})

add_executable(main ${SOURCE_FILES})
target_link_libraries(main m ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})