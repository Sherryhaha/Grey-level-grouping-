cmake_minimum_required(VERSION 3.6)
project(Grey_Level_grouping)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES  GLG.h gGLG.cpp GLG.cpp)
find_package(OpenCV REQUIRED)
add_executable(Grey_Level_grouping ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})