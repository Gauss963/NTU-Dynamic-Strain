cmake_minimum_required(VERSION 3.10)
project(Velocity-weakening LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TARGET_NAME Velocity-weakening)

find_package(Boost REQUIRED)
find_package(Akantu REQUIRED CONFIG)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)

# Akantu helper (adds proper include dirs/links for Akantu)
add_akantu_simulation(${TARGET_NAME} ${TARGET_NAME}.cc)

# Your extra deps
target_link_libraries(${TARGET_NAME} PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX)