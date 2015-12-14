find_package(Eigen3 REQUIRED)

include_directories(../exactdiag)

add_executable(free_fermion
               free_fermion.cc)

target_include_directories(free_fermion PRIVATE "${KORE_ROOT}" "${EIGEN3_INCLUDE_DIR}")
