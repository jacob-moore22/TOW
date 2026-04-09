cmake_minimum_required(VERSION 3.14)
project(sor_baseline CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(sor sor.cpp)
target_link_libraries(sor m)
