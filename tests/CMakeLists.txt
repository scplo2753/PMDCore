include(GoogleTest)

find_package(GTest REQUIRED)

set(PMDCORE_TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

set(PMDCORE_TEST_SOURCES
    LineParser_test.cpp
    Reconstruct_Alignment_test.cpp
    calPMD_test.cpp
    parsedData_class_tests.cpp
    statics_tests.cpp
    utilities_tests/utilities_math_test.cpp
    utilities_tests/utilities_test.cpp
)

add_executable(PMDCore_tests ${PMDCORE_TEST_SOURCES}
    statics_tests.cpp)

target_include_directories(PMDCore_tests
    PRIVATE
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/src"
)

target_link_libraries(PMDCore_tests
    PRIVATE
        SRC
        headers
        GTest::gtest_main
        Threads::Threads
        wrappers
)

target_compile_definitions(PMDCore_tests
    PRIVATE
        PMDCORE_TEST_DATA_DIR="${PMDCORE_TEST_DATA_DIR}"
)

gtest_discover_tests(PMDCore_tests
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    PROPERTIES
        ENVIRONMENT "PMDCORE_TEST_DATA_DIR=${PMDCORE_TEST_DATA_DIR}"
)
