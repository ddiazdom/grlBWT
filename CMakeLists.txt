cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
project(grlBWT VERSION 1.0.1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
find_package(LibSDSL REQUIRED)

# CMake options
option(MALLOC_COUNT "Compile grlBWT with malloc count (only for debugging purposes)" ON)
option(BUILD_SCRIPTS "Compile some helper scripts" ON)
option(MARCH_NATIVE "Enable -march=native when supported by the compiler" ON)
#

# load useful helper scripts
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils.cmake)
#

#get git info
get_git_commit(GIT_COMMIT)
#

#create header with program metadata
configure_file(
    ${CMAKE_SOURCE_DIR}/version.h.in
    ${CMAKE_BINARY_DIR}/generated/version.h
)
#

#make metadata header visible
include_directories(myprogram PRIVATE
    ${CMAKE_BINARY_DIR}/generated
)
#

add_executable(grlbwt-cli src/grlbwt-cli.cpp)
target_include_directories(grlbwt-cli PRIVATE include external)
configure_target(grlbwt-cli TRUE TRUE)

if (BUILD_SCRIPTS)
    function(create_script target_name source)
        add_executable(${target_name} ${source})
        target_include_directories(${target_name} PRIVATE include)
        configure_target(${target_name} TRUE TRUE)
    endfunction()
    create_script(reverse_bwt src/reverse_bwt.cpp)
    create_script(bwt_stats src/bwt_stats.cpp)
    create_script(split_runs src/split_runs.cpp)
    create_script(grl2rle src/grl2rle.cpp)
    create_script(grl2plain src/grl2plain.cpp)

    #testing
    create_script(test_new_rlbwt test_rlbwt.cpp)
endif ()