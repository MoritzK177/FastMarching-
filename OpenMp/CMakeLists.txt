cmake_minimum_required(VERSION 3.19)
project(FMMOpenMP)

set(CMAKE_CXX_STANDARD 20)


add_executable(FMMOpenMP main.cpp settings.h HeapAndStructVector.cpp HeapAndStructVector.h)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(FMMOpenMP PUBLIC OpenMP::OpenMP_CXX)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")