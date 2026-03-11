include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)

function(get_git_commit OUTPUT_VAR)
    find_package(Git QUIET)
    if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
        execute_process(
                COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_COMMIT
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        # Check if repo is dirty
        execute_process(
                COMMAND ${GIT_EXECUTABLE} diff --quiet
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                RESULT_VARIABLE GIT_DIRTY
        )
        if(NOT GIT_DIRTY EQUAL 0)
            set(GIT_COMMIT "${GIT_COMMIT}-dirty")
        endif()
    else()
        set(GIT_COMMIT "unknown")
    endif()
    set(${OUTPUT_VAR} "${GIT_COMMIT}" PARENT_SCOPE)
endfunction()

function(enable_march_native_if_requested target_name)
    if(MARCH_NATIVE)
        check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
        if(COMPILER_SUPPORTS_MARCH_NATIVE)
            message(STATUS "Compiler supports -march=native, enabling it for ${target_name}")
            target_compile_options(${target_name} PRIVATE -march=native)
        endif()
    else()
        message(STATUS "-march=native disabled for ${target_name}")
    endif()
endfunction()

function(enable_malloc_count_if_requested target_name)
    if(MALLOC_COUNT)
        message(STATUS "Compiling ${target_name} with malloc_count")
        target_sources(${target_name} PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR}/external/malloc_count-master/malloc_count.c
        )
        target_include_directories(${target_name} PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR}/external/malloc_count-master
        )
        target_compile_definitions(${target_name} PRIVATE USE_MALLOC_COUNT)

        if(UNIX AND NOT APPLE)
            target_link_libraries(${target_name} PRIVATE dl)
        endif()
    endif()
endfunction()

function(enable_filesystem_compat target_name)
    set(CMAKE_REQUIRED_QUIET TRUE)
    set(CMAKE_REQUIRED_LIBRARIES "")

    check_cxx_source_compiles("
        #include <filesystem>
        int main() {
            std::filesystem::path p{\".\"};
            return std::filesystem::exists(p) ? 0 : 0;
        }
    " FILESYSTEM_NO_EXTRA_LIB)

    if(NOT FILESYSTEM_NO_EXTRA_LIB)
        set(CMAKE_REQUIRED_LIBRARIES stdc++fs)

        check_cxx_source_compiles("
            #include <filesystem>
            int main() {
                std::filesystem::path p{\".\"};
                return std::filesystem::exists(p) ? 0 : 0;
            }
        " FILESYSTEM_WITH_STDCXXFS)

        if(FILESYSTEM_WITH_STDCXXFS)
            message(STATUS "std::filesystem requires linking stdc++fs for ${target_name}")
            target_link_libraries(${target_name} PRIVATE stdc++fs)
        else()
            message(FATAL_ERROR "std::filesystem is not usable with this compiler/toolchain")
        endif()
    endif()
endfunction()

function(configure_sdsl target_name)
    target_include_directories(${target_name} SYSTEM PRIVATE ${LIBSDSL_INCLUDE_DIRS})
    target_link_libraries(${target_name} PRIVATE ${LIBSDSL_LIBRARIES})
endfunction()

function(configure_target target_name enable_sdsl opt_flags)
    #standard flags
    target_compile_options(${target_name}
            PRIVATE
            -Wall -pedantic -Wno-deprecated-declarations
    )

    #optimization flags
    if(opt_flags)
        message(STATUS "Optimization flags enabled for ${target_name}")
        target_compile_options(${target_name}
                PRIVATE
                -O3 -funroll-loops -fomit-frame-pointer)
    else()
        message(STATUS "Optimization flags disabled for ${target_name}")
    endif()

    enable_march_native_if_requested(${target_name})
    enable_malloc_count_if_requested(${target_name})
    enable_filesystem_compat(${target_name})
    if(enable_sdsl)
        configure_sdsl(${target_name})
    endif()
endfunction()