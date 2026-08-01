file(GLOB SRC_FILES *.cpp)

add_library(SRC ${SRC_FILES}
    ThreadPool_entry.cpp
    statics.cpp)

target_include_directories(SRC PUBLIC ../include)

target_link_libraries(SRC PUBLIC Threads::Threads argparse)