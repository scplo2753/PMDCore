include(GoogleTest)

find_package(GTest REQUIRED)

set(PMDCORE_TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

add_executable(PMDCore_tests
    LineParser_test.cpp
    #Reconstruct_test.cpp
    parsedData_class_tests.cpp
)

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
