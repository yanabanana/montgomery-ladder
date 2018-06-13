cmake_minimum_required(VERSION 3.6)
project(montgomery)
include_directories(/usr/local/lib /usr/local/include)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp point.h)
find_library(MPIR_LIBRARY mpir)
find_library(MPIR_CPP_LIBRARY mpirxx)
add_executable(montgomery ${SOURCE_FILES})
target_link_libraries(montgomery "${MPIR_LIBRARY}" "${MPIR_CPP_LIBRARY}")
