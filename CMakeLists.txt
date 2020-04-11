cmake_minimum_required(VERSION 3.14)
project(DLAf_optimized)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(DLAf_optimized
        dlaf.cpp
        fastRandom.h
        nanoflann.hpp)
