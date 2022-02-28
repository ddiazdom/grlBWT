cmake_minimum_required(VERSION 3.7)

project(grl-BWT)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

find_package(LibSDSL REQUIRED)

add_executable(grlbwt
        main.cpp
        external/xxHash-dev/xxhash.c
        lib/grl_bwt.cpp
        lib/lc_gram_algo.cpp
        lib/utils.cpp
        lib/LMS_induction.cpp
        lib/fastx_handler.cpp
        lib/dna_string.cpp)

target_compile_options(grlbwt
        PRIVATE
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic
        -O3 -funroll-loops -fomit-frame-pointer -ffast-math
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

target_link_libraries(grlbwt LINK_PUBLIC z pthread ${LIBSDSL_LIBRARIES})

target_include_directories(grlbwt
        PRIVATE
        external
        ${LIBSDSL_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include)

add_executable(print_seqs scripts/print_seqs.cpp)

target_link_libraries(print_seqs
        LINK_PUBLIC
        pthread ${LIBSDSL_LIBRARIES})

target_include_directories(print_seqs
        PRIVATE
        ${LIBSDSL_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include)

add_executable(stats scripts/stats.cpp)
target_link_libraries(stats
        LINK_PUBLIC
        pthread ${LIBSDSL_LIBRARIES})

target_include_directories(stats
        PRIVATE
        ${LIBSDSL_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include)
#add_executable(tr
#        include/bwt_io.h
#        test_readers.cpp)
#
#target_compile_options(tr
#        PRIVATE
#        #-Werror
#        -Wall -Wextra -Wno-ignored-qualifiers -pedantic #-msse4.2 -O3 -funroll-loops -fomit-frame-pointer -ffast-math
#        INTERFACE
#        -Wshadow)
#
##no-stack-check only for APPLE!!!!!!
#if(APPLE)
#    target_compile_options(tr PRIVATE -fno-stack-check -Wno-vla-extension)
#endif()
#
#if(UNIX AND NOT APPLE)
#    target_link_libraries(tr LINK_PUBLIC stdc++fs)
#endif()
#
#target_link_libraries(tr LINK_PUBLIC pthread ${LIBSDSL_LIBRARIES})
#
#target_include_directories(tr
#        PRIVATE
#        external
#        ${LIBSDSL_INCLUDE_DIRS}
#        ${CMAKE_SOURCE_DIR}/include)
