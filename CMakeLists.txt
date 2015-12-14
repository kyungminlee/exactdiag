cmake_minimum_required(VERSION 3.3)
project(ed2)

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
        main.cc
)
set(HEADER_FILES
        exactdiag/operator.h
        exactdiag/hilbertspace.h exactdiag/quantumnumber.h)
include_directories(exactdiag)
add_executable(ed2 ${SOURCE_FILES} ${HEADER_FILES})


add_subdirectory(tests)