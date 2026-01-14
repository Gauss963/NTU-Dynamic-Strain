add_executable(frictionLaws
    main.cc
)

target_include_directories(frictionLaws
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# =========================
# ROOT (modern, clean)
# =========================
find_package(ROOT REQUIRED COMPONENTS
    Core
    Hist
    Graf
    Gpad
)

target_link_libraries(frictionLaws
    PRIVATE
        ROOT::Core
        ROOT::Hist
        ROOT::Graf
        ROOT::Gpad
)

target_compile_options(frictionLaws
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
)