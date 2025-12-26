cmake_minimum_required(VERSION 3.10)
project(Velocity-weakening-LSW LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TARGET_NAME Velocity-weakening-LSW)

find_package(Boost REQUIRED)
find_package(Akantu REQUIRED CONFIG)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)

add_akantu_simulation(${TARGET_NAME} ${TARGET_NAME}.cc)

target_link_libraries(${TARGET_NAME} PRIVATE
    /home/gauss112/softwares/scotch-v7.0.4/lib/libscotch.a
    /home/gauss112/softwares/scotch-v7.0.4/lib/libscotcherr.a
    /home/gauss112/softwares/scotch-v7.0.4/lib/libptscotch.a
    /home/gauss112/softwares/scotch-v7.0.4/lib/libptscotcherr.a
)