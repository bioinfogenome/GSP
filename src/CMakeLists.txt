cmake_minimum_required(VERSION 3.3)
project(GSP)

find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(GSP ${SOURCE_FILES})