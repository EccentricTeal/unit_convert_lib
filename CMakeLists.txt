cmake_minimum_required(VERSION 2.8.3)
project(unit_convert_lib)
add_compile_options(-std=c++14)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES unit_convert_lib
  # CATKIN_DEPENDS roscpp rospy
  # DEPENDS eigen
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)

add_library(unit_convert_lib SHARED
  src/coordinate.cc
  src/angle.cc
)
target_link_libraries(unit_convert_lib
  ${catkin_LIBRARIES}
)