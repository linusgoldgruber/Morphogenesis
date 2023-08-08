#cmake_minimum_required(VERSION 3.10)

## Set the project name and version
#project(DeviceInfo VERSION 1.0)

## Specify the C++ standard
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED True)



## Add the include directories
#include_directories(/usr/include/CL /home/goldi/Documents/KDevelop\ Projects/HandsOnOpenCL/Exercises-Solutions/Exercises/Cpp_common /home/goldi/Documents/KDevelop\ Projects/Morphogenesis/Morphogenesis/src)
#link_directories(/home/goldi/build/lib/libopencv_gapi.so.4.7.0)
##########################################################################

#find_package(OpenCL REQUIRED)
#find_package(jsoncpp REQUIRED)

##if( OpenCV_FOUND )
  ### Additional Include Directories
  ##include_directories( ${OpenCV_INCLUDE_DIRS} )

  ### Additional Library Directories
  ##link_directories( ${OpenCV_LIB_DIR} )

  ### Additional Dependencies

##endif()

###########################################################################

## Add an executable
#add_executable(DeviceInfo DeviceInfo.cpp)

#target_link_libraries(${PROJECT_NAME} PRIVATE
        #OpenCL::OpenCL
        #jsoncpp_lib
        #)

## Print the unclude directories used by the compiler
#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
      #foreach(dir ${dirs})
      #message(STATUS "dir='${dir}'")
      #endforeach()

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
cmake_policy(SET CMP0054 NEW)

project(DeviceInfo LANGUAGES C CXX)

set(CMAKE_C_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

find_package(OpenCL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

include_directories(${JSONCPP_INCLUDE_DIRS})

add_executable(DeviceInfo
DeviceInfo.cpp
)

target_link_libraries(${PROJECT_NAME} PRIVATE
        OpenCL::OpenCL
        ${JSONCPP_LIBRARIES}
)

set_target_properties(${PROJECT_NAME} PROPERTIES C_STANDARD 11
                                                 C_STANDARD_REQUIRED ON
                                                 C_EXTENSIONS OFF)

target_compile_definitions(${PROJECT_NAME} PRIVATE CL_TARGET_OPENCL_VERSION=300)
