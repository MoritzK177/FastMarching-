cmake_minimum_required(VERSION 3.19)
project(fastmarchingParallel1)

set(CMAKE_CXX_STANDARD 20)

add_executable(fastmarchingParallel1 main.cpp settings.h HeapAndStructVector.cpp HeapAndStructVector.h)