cmake_minimum_required(VERSION 3.14)

project(SplineConcepts)

# FetchContent for dependencies
include(FetchContent)

# NumCpp
FetchContent_Declare(
  numcpp
  GIT_REPOSITORY https://github.com/dpilger26/NumCpp.git
  GIT_TAG Version_2.14.0
)
FetchContent_MakeAvailable(numcpp)

# pybind11
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# tinyxml2
FetchContent_Declare(
  tinyxml2
  GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
  GIT_TAG 10.0.0
)
FetchContent_MakeAvailable(tinyxml2)

# Add library
file(GLOB_RECURSE SPLINECONCEPTS_SOURCES "${CMAKE_SOURCE_DIR}/splineconcepts/*.cpp" "${CMAKE_SOURCE_DIR}/splineconcepts/*.h")
# remove

# make export header
include(GenerateExportHeader)

add_library(splineconcepts ${SPLINECONCEPTS_SOURCES})
generate_export_header(splineconcepts)
target_include_directories(splineconcepts PUBLIC ${CMAKE_SOURCE_DIR}/splineconcepts ${NUMCPP_INCLUDE_DIRS} ${TINYXML2_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(splineconcepts PUBLIC NumCpp::NumCpp tinyxml2)

# Provide Python bindings
pybind11_add_module(pysplineconcepts ${CMAKE_SOURCE_DIR}/pysplineconcepts/python_bindings.cpp)
target_link_libraries(pysplineconcepts PRIVATE splineconcepts)

# Add tests
enable_testing()
file(GLOB_RECURSE TEST_SOURCES "${CMAKE_SOURCE_DIR}/tests/*.cpp")
foreach(TEST_SOURCE ${TEST_SOURCES})
  get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
  add_executable(${TEST_NAME} ${TEST_SOURCE})
  target_link_libraries(${TEST_NAME} splineconcepts)
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
