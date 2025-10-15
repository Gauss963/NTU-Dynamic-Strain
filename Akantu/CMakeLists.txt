cmake_minimum_required(VERSION 3.10)
project(cohesive_extrinsic_ig_tg)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(AKANTU_ROOT "/opt/homebrew/Caskroom/miniconda/base/envs/Akantu")

include_directories(${AKANTU_ROOT}/include)
include_directories(${AKANTU_ROOT}/include/akantu)
include_directories(${AKANTU_ROOT}/include/akantu_iterators)
include_directories(${AKANTU_ROOT}/include/eigen3)
include_directories(${AKANTU_ROOT}/include/iohelper)

# link_directories(${AKANTU_ROOT}/lib)
# add_executable(cohesive_extrinsic_ig_tg cohesive_extrinsic_ig_tg.cc)
# target_link_libraries(cohesive_extrinsic_ig_tg akantu iohelper)


add_executable(cohesive_extrinsic_ig_tg cohesive_extrinsic_ig_tg.cc)
target_link_libraries(cohesive_extrinsic_ig_tg
    ${AKANTU_ROOT}/lib/libakantu.dylib
    ${AKANTU_ROOT}/lib/libiohelper.dylib
)

target_include_directories(cohesive_extrinsic_ig_tg PRIVATE
    ${AKANTU_ROOT}/include
    ${AKANTU_ROOT}/include/akantu
    ${AKANTU_ROOT}/include/akantu_iterators
    ${AKANTU_ROOT}/include/eigen3
    ${AKANTU_ROOT}/include/iohelper
)

set_target_properties(cohesive_extrinsic_ig_tg PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)



target_compile_features(cohesive_extrinsic_ig_tg PRIVATE cxx_std_17)