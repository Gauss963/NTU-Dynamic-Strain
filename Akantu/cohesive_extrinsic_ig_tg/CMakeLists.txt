cmake_minimum_required(VERSION 3.10)
# project(cohesive_extrinsic_ig_tg)
project(Velocity-weakening)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Akantu path (your conda env)
set(AKANTU_ROOT "/opt/homebrew/Caskroom/miniconda/base/envs/Akantu")
# set(TARGET_NAME cohesive_extrinsic_ig_tg)
set(TARGET_NAME Velocity-weakening)

# ───── Include Akantu headers ─────
include_directories(${AKANTU_ROOT}/include)
include_directories(${AKANTU_ROOT}/include/akantu)
include_directories(${AKANTU_ROOT}/include/akantu/materials)
include_directories(${AKANTU_ROOT}/include/akantu_iterators)
include_directories(${AKANTU_ROOT}/include/eigen3)
include_directories(${AKANTU_ROOT}/include/iohelper)

# ───── Link Akantu libraries ─────
link_directories(${AKANTU_ROOT}/lib)

# ───── Find OpenMP ─────
find_package(OpenMP REQUIRED COMPONENTS CXX)

# ───── Build target ─────
add_executable(${TARGET_NAME} ${TARGET_NAME}.cc)

# ───── Link libraries ─────
target_link_libraries(${TARGET_NAME}
    akantu
    iohelper
    OpenMP::OpenMP_CXX
)

# ───── Compile features ─────
target_compile_features(${TARGET_NAME} PRIVATE cxx_std_17)