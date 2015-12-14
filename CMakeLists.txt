cmake_minimum_required(VERSION 3.3)
project(ed2)


get_filename_component(KORE_ROOT "../../../koretools/include" ABSOLUTE)

set(KORE_ROOT ../../../koretools/include)

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

set(SOURCE_FILES
    main.cc
    )
set(HEADER_FILES
    exactdiag/operator.h
    exactdiag/hilbertspace.h exactdiag/hilbertspace/quantumnumber.h exactdiag/operator/pure_operator.h exactdiag/operator/raw_rep_operator.h examples/free_fermion.cc exactdiag/hilbertspace/state.h exactdiag/hilbertspace/site.h exactdiag/hilbertspace/system.h exactdiag/hilbertspace/basis_iterator.h exactdiag/hilbertspace/sector.h exactdiag/global.h exactdiag/operator/generic_operator.h)

include_directories(exactdiag "${KORE_ROOT}")
add_executable(ed2 ${SOURCE_FILES} ${HEADER_FILES})


add_subdirectory(tests)
add_subdirectory(examples)