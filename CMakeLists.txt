cmake_minimum_required(VERSION 3.15)
project(monod_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Enable testing
enable_testing()

# Build the library
add_library(monod_model STATIC monod_model.cpp)
target_include_directories(monod_model PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(monod_model PUBLIC Eigen3::Eigen)

# Add validation tests
add_executable(test_monod_validation test_monod_cpp_main.cpp)
target_link_libraries(test_monod_validation PRIVATE monod_model)

add_executable(test_reactive_transport_validation test_reactive_transport_2d.cpp)
target_link_libraries(test_reactive_transport_validation PRIVATE monod_model Eigen3::Eigen)

# Register tests with CTest
add_test(
    NAME monod_validation
    COMMAND test_monod_validation
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME reactive_transport_validation
    COMMAND test_reactive_transport_validation
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Set test properties
set_tests_properties(
    monod_validation
    reactive_transport_validation
    PROPERTIES
    TIMEOUT 30
    FAIL_REGULAR_EXPRESSION "FAILED|Error|error"
)
