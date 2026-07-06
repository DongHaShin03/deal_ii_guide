cmake_minimum_required(VERSION 3.12)

project(step-1 LANGUAGES CXX)

include(${CMAKE_CURRENT_SOURCE_DIR}/../1_common/cmake-common.cmake)

add_executable(step_1 src/step_1.cpp)

deal_ii_setup_target(step_1)
