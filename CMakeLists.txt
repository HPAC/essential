cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(MC-GEN LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(comp_flags -Wall)

# Enable different parts
option(WITH_TESTS "Enable Tests" YES)

add_subdirectory(src)

if(WITH_TESTS)
  add_subdirectory(test)
endif()