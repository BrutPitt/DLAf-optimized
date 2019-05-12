# CMakeList.txt : CMake project for dlaf, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (dlaf "../dlaf.cpp" "../nanoflann.hpp" "../fastRandom.h")

# TODO: Add tests and install targets if needed.
