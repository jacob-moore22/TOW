cmake_minimum_required(VERSION 3.18)
project(sor_matar_serial CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

FetchContent_Declare(
    matar
    GIT_REPOSITORY https://github.com/lanl/MATAR.git
    GIT_TAG        main
    GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(matar)
if(NOT matar_POPULATED)
    FetchContent_Populate(matar)
endif()

add_library(matar_headers INTERFACE)
target_include_directories(matar_headers INTERFACE ${matar_SOURCE_DIR}/src/include)

add_executable(sor sor.cpp)
target_link_libraries(sor matar_headers)
