cmake_minimum_required(VERSION 3.7)

project(grlBWT)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

find_package(LibSDSL REQUIRED)
add_subdirectory(external/cdt)
add_subdirectory(external/bioparsers)

add_library(grlbwt
        lib/grl_bwt.cpp
        lib/lc_gram_algo.cpp
        lib/LMS_induction.cpp)

target_compile_options(grlbwt
        PRIVATE
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic
        #-O3 -funroll-loops -fomit-frame-pointer -ffast-math
        INTERFACE
        -Wshadow)

# this extra flags only applies to clang compilers
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Disabling vla-extension due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(grlbwt PUBLIC -Wno-vla-extension -Wno-undefined-var-template)
endif()

# there seems to be problem with the msse4.2 compiler flag and the new Apple chips
if(NOT CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64")
    target_compile_options(grlbwt PUBLIC -msse4.2)
else()
    message(STATUS "Disabling SSE4.2 instructions due to conflict with host system processor ${CMAKE_HOST_SYSTEM_PROCESSOR}")
endif()

#for the file system functionality on Linux
if(UNIX AND NOT APPLE)
    target_link_libraries(grlbwt LINK_PUBLIC stdc++fs)
endif()

target_link_libraries(grlbwt LINK_PUBLIC z pthread ${LIBSDSL_LIBRARIES} cdt)

target_include_directories(grlbwt
        PUBLIC
        ${LIBSDSL_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/external)

target_include_directories(grlbwt SYSTEM PUBLIC ${LIBSDSL_INCLUDE_DIRS})

add_executable(grlbwt-cli main.cpp)
target_link_libraries(grlbwt-cli grlbwt biopar)
target_include_directories(grlbwt-cli PUBLIC external/bioparsers/include)

