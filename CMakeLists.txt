cmake_minimum_required(VERSION 3.8)
project(strategy-optimizer)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES strategy.cpp)
add_executable(strategy-optimizer ${SOURCE_FILES})