cmake_minimum_required(VERSION 3.5)
project(unit_convert_lib)

# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(Eigen3 REQUIRED)



# add_library(unit_convert_lib
#   INTERFACE #HEADER ONLY
# )
add_library(convert_angle SHARED
  src/angle.cc
)
target_compile_definitions(convert_angle
  PRIVATE "CONVERT_ANGLE_BUILDING_LIBRARY"
)
ament_target_dependencies(convert_angle
  Eigen3
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
ament_target_dependencies(convert_coordinate
  Eigen3
)
target_include_directories(convert_coordinate
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${Eigen3_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)



if(NOT WIN32)
  ament_environment_hooks(
    "${ament_cmake_package_templates_ENVIRONMENT_HOOK_LIBRARY_PATH}"
  )
endif()
ament_export_targets(export_${PROJECT_NAME} HAS_LIBRARY_TARGET)

install(DIRECTORY include/
  DESTINATION include
)
install(TARGETS convert_angle convert_coordinate
  EXPORT export_${PROJECT_NAME}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)


ament_package()
