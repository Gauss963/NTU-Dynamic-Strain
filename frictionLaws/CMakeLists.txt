cmake_minimum_required(VERSION 3.16)

project(frictionLaws
    VERSION 0.1
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ---- ROOT auto-detection safeguard ----
if(NOT ROOT_DIR AND DEFINED ENV{ROOTSYS})
    set(ROOT_DIR "$ENV{ROOTSYS}/cmake")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_subdirectory(src)