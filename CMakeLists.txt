cmake_minimum_required(VERSION 3.16.3)
project(Raytracing)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(glm)

add_executable(Raytracing
        main.cpp
        geometry/geometry.cpp
        geometry/geometry.h
        geometry/scene.h
        distributions/distributions.cpp
        distributions/distributions.h
        distributions/engine.h)

target_include_directories(Raytracing PUBLIC "${CMAKE_CURRENT_LIST_DIR}")

target_link_libraries(Raytracing PUBLIC glm)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
    message(FATAL_ERROR "OpenMP not found! You need OpenMP for speedup on multicore CPUs!")
endif()
