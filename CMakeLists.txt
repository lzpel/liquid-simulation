cmake_minimum_required(VERSION 3.12)
project(main)

set(MIN_COMPILE "-O2 -Os -DNDEBUG -s")

include_directories(src/)
file(GLOB_RECURSE src "src/*.cpp" "src/*.c")
add_executable(${PROJECT_NAME} ${SYSTEM} ${src} ${dep})
target_link_libraries(${PROJECT_NAME})