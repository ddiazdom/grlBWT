cmake_minimum_required(VERSION 3.7)

project(grlBWT-cli)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)
add_executable(grlbwt-cli main.cpp)
target_link_libraries(grlbwt-cli grlbwt)

#
add_executable(reverse_bwt scripts/reverse_bwt.cpp)
target_link_libraries(reverse_bwt grlbwt)
target_include_directories(reverse_bwt PUBLIC scripts/)
#
#add_executable(bwt_stats scripts/bwt_stats.cpp)
#target_link_libraries(bwt_stats grlbwt biopar)
#target_include_directories(bwt_stats PUBLIC scripts/)
#
#add_executable(split_runs scripts/split_runs.cpp)
#target_link_libraries(split_runs grlbwt biopar)
#target_include_directories(split_runs PUBLIC scripts/)
#
#add_executable(grlbwt2rle scripts/grlbwt2rle.cpp)
#target_link_libraries(grlbwt2rle grlbwt biopar)
#target_include_directories(grlbwt2rle PUBLIC scripts/)
#
#add_executable(grl2plain scripts/grl2plain.cpp)
#target_link_libraries(grl2plain grlbwt biopar)
#target_include_directories(grl2plain PUBLIC scripts/)