cmake_minimum_required(VERSION 3.15)
project(PMDCore VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Only enable tests for Debug builds.
include(CTest)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(__DEBUG__)
    set(BUILD_TESTING ON CACHE BOOL "Enable tests" FORCE)
else()
    set(BUILD_TESTING OFF CACHE BOOL "Enable tests" FORCE)
endif()
#set(CMAKE_UNITY_BUILD ON)

find_package(Threads REQUIRED)

add_executable(PMDCore main_multithreaded.cpp)

option(PMDCORE_STATIC_DEPS "Prefer static third-party/internal libraries" OFF)

if(PMDCORE_STATIC_DEPS)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build libraries as shared objects" FORCE)
else()
    set(BUILD_SHARED_LIBS ON CACHE BOOL "Build libraries as shared objects" FORCE)
endif()

add_subdirectory(third_party/matplotplusplus)
add_subdirectory(third_party/argparse)
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(wrappers)

target_link_libraries(PMDCore PRIVATE SRC 
matplot argparse Threads::Threads wrappers)

if(PMDCORE_STATIC_DEPS AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_options(PMDCore PRIVATE -static-libgcc -static-libstdc++)
endif()

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
