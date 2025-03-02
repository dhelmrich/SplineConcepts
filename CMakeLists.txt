cmake_minimum_required(VERSION 3.14)

project(SplineConcepts)

# Set default build folder if not otherwise specified
if(NOT CMAKE_BUILD_TYPE)
  if(WIN32)
    set(CMAKE_BUILD_TYPE "build_windows")
  elseif(UNIX AND NOT APPLE)
    set(CMAKE_BUILD_TYPE "build_linux")
  elseif(APPLE)
    set(CMAKE_BUILD_TYPE "build_macos")
  endif()
endif()


# FetchContent for dependencies
include(FetchContent)

# set CMP0135 policy
cmake_policy(SET CMP0135 NEW)

# Set to use C++17
set(CMAKE_CXX_STANDARD 17)

# NumCpp
# Declare NumCPP options:
option(NUMCPP_NO_USE_BOOST "Do not use boost" ON)
option(NUMCPP_USE_MULTITHREAD "Use Multithreding" ON)

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

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

# Add library
file(GLOB_RECURSE SPLINECONCEPTS_SOURCES "${CMAKE_SOURCE_DIR}/splineconcepts/*.cpp" "${CMAKE_SOURCE_DIR}/splineconcepts/*.hpp")
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
  target_link_libraries(${TEST_NAME} splineconcepts GTest::gtest_main)
  target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/tests ${NUMCPP_INCLUDE_DIRS} ${TINYXML2_INCLUDE_DIRS})
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
  gtest_discover_tests(${TEST_NAME})
endforeach()

# Add examples
file(GLOB_RECURSE EXAMPLE_SOURCES "${CMAKE_SOURCE_DIR}/examples/*.cpp")
foreach(EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
  get_filename_component(EXAMPLE_NAME ${EXAMPLE_SOURCE} NAME_WE)
  add_executable(${EXAMPLE_NAME} ${EXAMPLE_SOURCE})
  target_link_libraries(${EXAMPLE_NAME} splineconcepts)
  target_include_directories(${EXAMPLE_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/examples ${NUMCPP_INCLUDE_DIRS} ${TINYXML2_INCLUDE_DIRS})
endforeach()

