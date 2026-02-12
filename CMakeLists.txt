cmake_minimum_required(VERSION 3.15)
project(monod_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(monod_model STATIC monod_model.cpp)
target_include_directories(monod_model PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_monod_cpp_main test_monod_cpp_main.cpp)
target_link_libraries(test_monod_cpp_main PRIVATE monod_model)
