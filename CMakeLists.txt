cmake_minimum_required(VERSION 3.25)

project(CRTFactoring LANGUAGES CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(GMP_ROOT "C:/Strawberry/c")
find_path(GMP_INCLUDE_DIR gmp.h PATHS "${GMP_ROOT}/include" NO_DEFAULT_PATH)
find_file(GMP_DLL NAMES libgmp-10.dll PATHS "${GMP_ROOT}/bin" NO_DEFAULT_PATH)

if(NOT GMP_INCLUDE_DIR OR NOT GMP_DLL)
    message(FATAL_ERROR "GMP was not found under ${GMP_ROOT}. Install GMP or update GMP_ROOT in CMakeLists.txt.")
endif()

set(GMP_IMPORT_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/gmp-msvc.lib")
set(GMP_EXPORTS_FILE "${CMAKE_CURRENT_SOURCE_DIR}/gmp/gmp-msvc.def")
add_custom_command(
    OUTPUT "${GMP_IMPORT_LIBRARY}"
    COMMAND "${CMAKE_AR}"
            "/def:${GMP_EXPORTS_FILE}"
            /machine:x64
            "/out:${GMP_IMPORT_LIBRARY}"
    DEPENDS "${GMP_EXPORTS_FILE}"
    COMMENT "Creating the MSVC import library for GMP"
    VERBATIM
)
add_custom_target(gmp_import_library DEPENDS "${GMP_IMPORT_LIBRARY}")

function(add_factor_app target)
    add_executable(${target} ${ARGN})
    add_dependencies(${target} gmp_import_library)
    target_include_directories(${target} PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/gmp"
        "${GMP_INCLUDE_DIR}"
    )
    target_link_libraries(${target} PRIVATE "${GMP_IMPORT_LIBRARY}")
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${GMP_DLL}"
                "$<TARGET_FILE_DIR:${target}>"
        VERBATIM
    )
endfunction()

add_factor_app(factorbase
    src/factorbase.cpp
    src/math.cpp
    src/crt.cpp
    src/main_crt.cpp
)

# Optimized executable intended for ordinary experiment runs. Keep this as a
# separate target so Debug builds retain their useful checks and symbols.
add_factor_app(factorbase_optimized
    src/factorbase.cpp
    src/math.cpp
    src/crt.cpp
    src/main_crt.cpp
)
if(MSVC)
    target_compile_options(factorbase_optimized PRIVATE /O2 /Ob2 /Oi /Ot /GL /DNDEBUG)
    target_link_options(factorbase_optimized PRIVATE /LTCG)
else()
    target_compile_options(factorbase_optimized PRIVATE -O3 -DNDEBUG)
endif()

add_factor_app(cofactor_experiment
    src/math.cpp
    experiments/fixed_cofactor.cpp
)

add_factor_app(cofactor_intersections
    src/math.cpp
    experiments/cofactor_intersections.cpp
)
find_package(OpenMP REQUIRED COMPONENTS CXX)
target_link_libraries(cofactor_intersections PRIVATE OpenMP::OpenMP_CXX)

if(BUILD_TESTING)
    add_factor_app(crt_test
        src/factorbase.cpp
        src/math.cpp
        src/crt.cpp
        tests/crt_test.cpp
    )
    add_test(NAME crt_refinement_vectors COMMAND crt_test)
    add_test(NAME fixed_cofactor_invariants COMMAND cofactor_experiment --self-test)
    add_test(NAME cofactor_intersection_invariants COMMAND cofactor_intersections --self-test)
    add_test(NAME cofactor_openmp_reproducibility COMMAND cofactor_intersections --parallel-self-test)
endif()
