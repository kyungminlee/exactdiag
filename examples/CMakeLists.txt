
include_directories(../exactdiag)

add_executable(free_fermion
               free_fermion.cc)

target_include_directories(free_fermion PRIVATE "${KORE_ROOT}" "eigen")

