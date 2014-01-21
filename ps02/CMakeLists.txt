project(ps02)
enable_language(CXX)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Weffc++")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O2")

set(SRC code.cpp)

add_executable(code ${SRC})