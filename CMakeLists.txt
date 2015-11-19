cmake_minimum_required(VERSION 3.3)
project(ed2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
        operator.h
        operator_test.cc hilbertspace.h quantumnumber.h)

add_executable(ed2 ${SOURCE_FILES})