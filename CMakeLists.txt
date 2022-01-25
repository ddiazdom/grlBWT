cmake_minimum_required(VERSION 3.7)
project(gBWT)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

find_package(LibSDSL REQUIRED)

add_executable(gbwt
        main.cpp
        external/xxHash-dev/xxhash.c
        lib/grammar.cpp
        lib/gbwt.cpp
        lib/lc_gram_algo.cpp
        lib/utils.cpp
        lib/LMS_induction.cpp)

target_compile_options(gbwt
        PRIVATE
        #-Werror
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic #-msse4.2 -O3 -funroll-loops -fomit-frame-pointer -ffast-math
        INTERFACE
        -Wshadow)

#no-stack-check only for APPLE!!!!!!
if(APPLE)
    target_compile_options(gbwt PRIVATE -fno-stack-check -Wno-vla-extension)
endif()

if(UNIX AND NOT APPLE)
    target_link_libraries(gbwt LINK_PUBLIC stdc++fs)
endif()

target_link_libraries(gbwt LINK_PUBLIC pthread ${LIBSDSL_LIBRARIES})

target_include_directories(gbwt
        PRIVATE
        external
        ${LIBSDSL_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include)

#add_executable(tr
#        include/bwt_io.h
#        test_readers.cpp)

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
