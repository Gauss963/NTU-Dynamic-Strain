cmake_minimum_required(VERSION 3.10)

project(Velocity-weakening-LSW LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TARGET_NAME Velocity-weakening-LSW)

# ---- Dependencies ----
find_package(Boost REQUIRED)
find_package(Akantu REQUIRED CONFIG)

# OpenMP is optional but allowed
find_package(OpenMP QUIET COMPONENTS CXX)

# ---- Akantu simulation target ----
add_akantu_simulation(${TARGET_NAME} ${TARGET_NAME}.cc)

# ---- OpenMP (only if available) ----
if(OpenMP_CXX_FOUND)
  target_link_libraries(${TARGET_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()