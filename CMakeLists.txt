cmake_minimum_required(VERSION 3.7)
project(LPG_COMPRESSOR)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

find_package(LibSDSL REQUIRED)

add_executable(lpg
        main.cpp
        external/xxHash-dev/xxhash.c
        lib/grammar.cpp
        lib/grammar_build.cpp
        lib/lc_gram_algo.cpp
        lib/utils.cpp
        lib/suffpair_algo.cpp include/lc_parsings.hpp)

target_compile_options(lpg
        PRIVATE
        #-Werror
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic -Wno-vla-extension -msse4.2 -O3 -funroll-loops -fomit-frame-pointer -ffast-math
        INTERFACE
        -Wshadow)

#no-stack-check only for APPLE!!!!!!
if(APPLE)
    target_compile_options(lpg PRIVATE -fno-stack-check)
endif()

target_link_libraries(lpg
        LINK_PUBLIC
        pthread ${LIBSDSL_LIBRARIES})

target_include_directories(lpg
        PRIVATE
        external
        ${LIBSDSL_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include)

