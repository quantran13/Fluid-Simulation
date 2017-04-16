cmake_minimum_required(VERSION 3.5)
project(fluid_simulation)

#set(CMAKE_C_STANDARD 11)

set(SOURCE_FILES src/main.cpp src/graphic.cu src/fluid.cu includes/fluid.h
        includes/graphic.h src/utility.cu includes/utility.h
        includes/fluid_kernels.h src/fluid_kernels.cu)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(CUDA QUIET REQUIRED)

include_directories(includes/)
include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})

cuda_add_executable(main ${SOURCE_FILES})
target_link_libraries(main m ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})