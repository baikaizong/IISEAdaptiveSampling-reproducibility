cmake_minimum_required(VERSION 3.25)

project(mMKSTD_reproducibility VERSION 1.0.0 LANGUAGES Fortran CXX)

set(IMSL_ROOT "$ENV{IMSL_ROOT}" CACHE PATH
    "Root of the architecture-specific IMSL Fortran Numerical Library installation")
set(GUROBI_ROOT "$ENV{GUROBI_ROOT}" CACHE PATH "Root of the Gurobi installation")

if(NOT GUROBI_ROOT AND DEFINED ENV{GUROBI_HOME})
    set(GUROBI_ROOT "$ENV{GUROBI_HOME}" CACHE PATH "Root of the Gurobi installation" FORCE)
endif()

if(NOT IMSL_ROOT)
    message(FATAL_ERROR "Set IMSL_ROOT to the architecture-specific IMSL installation directory.")
endif()
if(NOT GUROBI_ROOT)
    message(FATAL_ERROR "Set GUROBI_ROOT to the Gurobi installation directory.")
endif()

find_package(MPI REQUIRED COMPONENTS Fortran)

find_path(IMSL_INCLUDE_DIR RNSET_INT.mod
    HINTS "${IMSL_ROOT}/include/dll" "${IMSL_ROOT}/include/static"
    NO_DEFAULT_PATH)
find_library(IMSL_LIBRARY NAMES imsl_dll imsl
    HINTS "${IMSL_ROOT}/lib"
    NO_DEFAULT_PATH)
if(NOT IMSL_INCLUDE_DIR OR NOT IMSL_LIBRARY)
    message(FATAL_ERROR "IMSL was not found. Set IMSL_ROOT to the architecture-specific IMSL directory.")
endif()

add_library(project_dependencies INTERFACE)
target_include_directories(project_dependencies INTERFACE "${IMSL_INCLUDE_DIR}")
target_link_libraries(project_dependencies INTERFACE
    MPI::MPI_Fortran
    "${IMSL_LIBRARY}")

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(project_dependencies INTERFACE
        $<$<COMPILE_LANGUAGE:Fortran>:/warn:interfaces>
        $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:/Od>
        $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:/traceback>
        $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:/check:bounds>
        $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:/check:stack>
        $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:/heap-arrays:0>)
endif()

find_path(GUROBI_INCLUDE_DIR gurobi_c.h HINTS "${GUROBI_ROOT}/include" NO_DEFAULT_PATH)
find_library(GUROBI_LIBRARY NAMES gurobi130 gurobi HINTS "${GUROBI_ROOT}/lib" NO_DEFAULT_PATH)
if(NOT GUROBI_INCLUDE_DIR OR NOT GUROBI_LIBRARY)
    message(FATAL_ERROR "Gurobi was not found under GUROBI_ROOT.")
endif()

add_library(gurobi_wrapper STATIC src/gurobi/gurobi_interface.cpp)
target_include_directories(gurobi_wrapper PRIVATE "${GUROBI_INCLUDE_DIR}")
target_link_libraries(gurobi_wrapper PUBLIC "${GUROBI_LIBRARY}")

set(CORE_SOURCES
    src/core/GlobalSettings_mod.f90
    src/core/SparseMatrix_mod.f90
    src/core/Bsplines_mod.f90
    src/core/utils_mod.f90
    src/core/SamplesGeneration_mod.f90
    src/core/mMSTD_mod.f90
    src/core/CDS_mod.f90
    src/core/NAS_mod.f90
    src/core/POS_mod.f90
    src/core/SASAM_mod.f90
    src/core/TOPR_mod.f90
    src/core/TSBSS_mod.f90
    src/core/TSSPR_mod.f90
    src/core/Exploration_mod.f90
    src/core/Data_Loader.f90
    src/core/Gurobi_Interface_mod.f90)

add_library(methods STATIC ${CORE_SOURCES})
set(methods_module_dir "${CMAKE_CURRENT_BINARY_DIR}/modules/methods")
set_target_properties(methods PROPERTIES Fortran_MODULE_DIRECTORY "${methods_module_dir}")
target_include_directories(methods PUBLIC "$<BUILD_INTERFACE:${methods_module_dir}>")
target_link_libraries(methods PUBLIC project_dependencies gurobi_wrapper)

function(add_evaluation target app_dir)
    add_executable(${target}
        src/support/RuntimeConfig_mod.f90
        src/support/ExperimentSupport_mod.f90
        "src/apps/${app_dir}/performance_evaluation.f90"
        "src/apps/${app_dir}/Experiment_Runner_mod.f90"
        "src/apps/${app_dir}/main.f90")
    set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/modules/${target}")
    set_target_properties(${target} PROPERTIES
        Fortran_MODULE_DIRECTORY "${module_dir}"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
    target_include_directories(${target} PRIVATE "${module_dir}")
    target_link_libraries(${target} PRIVATE methods project_dependencies)
endfunction()

add_executable(case_study CaseStudy/CaseStudy.f90)
set(case_study_module_dir "${CMAKE_CURRENT_BINARY_DIR}/modules/case_study")
set_target_properties(case_study PROPERTIES
    Fortran_MODULE_DIRECTORY "${case_study_module_dir}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
target_include_directories(case_study PRIVATE "${case_study_module_dir}")
target_link_libraries(case_study PRIVATE project_dependencies)
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(case_study PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:/heap-arrays:0>)
endif()

add_evaluation(mmkstd_evaluation mmkstd)
add_evaluation(cds_evaluation cds)
add_evaluation(nas_evaluation nas)
add_evaluation(pos_evaluation pos)
add_evaluation(sasam_evaluation sasam)
add_evaluation(topr_evaluation topr)
add_evaluation(tsbss_evaluation tsbss)
add_evaluation(tsspr_evaluation tsspr)
