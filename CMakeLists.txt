
cmake_minimum_required(VERSION 3.12)

project(hypercutter LANGUAGES CXX VERSION 0.0.1)


set (CMAKE_CXX_STANDARD 20)
set (MYAPP leetcode)



include_directories(third_party/eigen-3.4.0/)
include_directories(third_party/slimcpplib/include)
include_directories(src/bsp/)


add_subdirectory(main)

# message(STATUS "Looking for someheader.h - found" ${CMAKE_SOURCE_DIR})

