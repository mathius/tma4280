project(hello)
enable_language(CXX)
cmake_minimum_required(VERSION 2.6)

set(SRC matrix.cpp)

add_executable(matrix ${SRC})