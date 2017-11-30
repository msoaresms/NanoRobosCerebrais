cmake_minimum_required(VERSION 3.8)
project(dontTalkToMeAboutLife)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp)
add_executable(dontTalkToMeAboutLife ${SOURCE_FILES})

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")