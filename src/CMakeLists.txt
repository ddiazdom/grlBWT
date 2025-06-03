project(grlBWT)

set(CMAKE_CXX_STANDARD 17)

option(MALLOC_COUNT "Compiling grlBWT with malloc count (only for debugging purposes)" OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

find_package(LibSDSL REQUIRED)

set(GRLBWT_LIB_FILES
        lib/ind_phase.cpp
        lib/par_phase.cpp
        lib/xxhash.c
        lib/utils.cpp
)

if(MALLOC_COUNT)
    list(APPEND GRLBWT_LIB_FILES malloc_count-master/malloc_count.c)
endif()

add_library(grlbwt ${GRLBWT_LIB_FILES})

if(MALLOC_COUNT)
    target_compile_definitions(grlbwt PRIVATE USE_MALLOC_COUNT)
endif()

target_compile_options(grlbwt
        PRIVATE
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic
        -O3 -funroll-loops -fomit-frame-pointer -ffast-math
        INTERFACE
        -Wshadow
)

#compiler-dependent flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Disabling vla-extension due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(grlbwt PUBLIC -Wno-vla-extension -Wno-undefined-var-template)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(grlbwt PUBLIC -march=native)
endif()

if(NOT CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64")
    target_compile_options(grlbwt PUBLIC -msse4.2)
endif()

#for the file_system functionality on Linux
if(UNIX AND NOT APPLE)
    target_link_libraries(grlbwt LINK_PUBLIC stdc++fs)
endif()

if(MALLOC_COUNT)
    message(STATUS "Adding malloc count")
    target_link_libraries(grlbwt LINK_PUBLIC dl)
    target_include_directories(grlbwt SYSTEM PUBLIC malloc_count-master)
endif()

target_link_libraries(grlbwt LINK_PUBLIC pthread z ${LIBSDSL_LIBRARIES})
target_include_directories(grlbwt PUBLIC include)
target_include_directories(grlbwt SYSTEM PUBLIC ${LIBSDSL_INCLUDE_DIRS})