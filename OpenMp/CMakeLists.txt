cmake_minimum_required(VERSION 3.19)
project(FMMOpenMP)

set(CMAKE_CXX_STANDARD 20)


add_executable(FMMOpenMP main.cpp settings.h HeapAndStructVector.cpp HeapAndStructVector.h)

