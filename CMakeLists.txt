cmake_minimum_required(VERSION 2.8.3)
project(unit_convert_lib)

find_package(catkin REQUIRED COMPONENTS
  Eigen
)

catkin_python_setup()

catkin_package(
  INCLUDE_DIRS include
#  LIBRARIES shipcon
#  CATKIN_DEPENDS roscpp rospy rqt_gui rqt_gui_py std_msgs message_runtime
#  DEPENDS system_lib
)

###########
## Build ##
###########
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)


# add_library(unit_convert_lib
#   INTERFACE #HEADER ONLY
# )
add_library(convert_angle SHARED
  src/angle.cc
)
target_compile_definitions(convert_angle
  PRIVATE "CONVERT_ANGLE_BUILDING_LIBRARY"
)
target_include_directories(convert_angle
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)



add_library(convert_coordinate SHARED
  src/coordinate.cc
)
target_compile_definitions(convert_coordinate
  PRIVATE "CONVERT_COORDINATE_BUILDING_LIBRARY"
)
target_include_directories(convert_coordinate
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${Eigen3_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)


#############
## Install ##
#############
 install(FILES
   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
 )
