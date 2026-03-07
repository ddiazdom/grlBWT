cmake_minimum_required(VERSION 3.10)

project(grlBWT)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
find_package(LibSDSL REQUIRED)

# CMake options
option(MALLOC_COUNT "Compile grlBWT with malloc count (only for debugging purposes)" OFF)
option(BUILD_SCRIPTS "Compile some helper scripts" ON)
option(MARCH_NATIVE "Enable -march=native when supported by the compiler" ON)
#

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils.cmake)

add_executable(grlbwt-cli main.cpp)
target_include_directories(grlbwt-cli PRIVATE include external)
configure_target(grlbwt-cli TRUE TRUE)

if (BUILD_SCRIPTS)
    function(create_script target_name source)
        add_executable(${target_name} ${source})
        target_include_directories(${target_name} PRIVATE scripts include)
        configure_target(${target_name} TRUE TRUE)
    endfunction()
    create_script(reverse_bwt scripts/reverse_bwt.cpp)
    create_script(bwt_stats scripts/bwt_stats.cpp)
    create_script(split_runs scripts/split_runs.cpp)
    create_script(grl2rle scripts/grl2rle.cpp)
    create_script(grl2plain scripts/grl2plain.cpp)
endif ()