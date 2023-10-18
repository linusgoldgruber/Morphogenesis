cmake_minimum_required(VERSION 3.10)
set(GCC_COVERAGE_COMPILE_FLAGS "-m32")
cmake_policy(SET CMP0054 NEW)

project(MinimalSim LANGUAGES C CXX)

set(PROJECT_NAME check_demo)
set(CMAKE_C_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I /usr/lib/gcc/x86_64-linux-gnu/11/include -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
set(VTK_DIR "/home/goldi/Documents/VTK/build")
set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

find_package(Boost  REQUIRED COMPONENTS filesystem) # see https://cliutils.gitlab.io/modern-cmake/chapters/packages/Boost.html
message(STATUS "Boost version: ${Boost_VERSION}")

# This is needed if your Boost version is newer than your CMake version
# or if you have an old version of CMake (<3.5)
if(NOT TARGET Boost::filesystem)
    add_library(Boost::filesystem IMPORTED INTERFACE)
    set_property(TARGET Boost::filesystem PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::filesystem PROPERTY
        INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

find_package(OpenCL REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(jsoncpp REQUIRED)
find_package(Threads REQUIRED)
find_package(VTK REQUIRED)

pkg_check_modules(JSONCPP jsoncpp)

include_directories(${JSONCPP_INCLUDE_DIRS})
include_directories(${VTK_INCLUDE_DIRS})

add_executable(${PROJECT_NAME}
        check_demo.cpp
        )

#add_executable(CreateDemo
        #createDemo.cpp
        #)

target_link_libraries(${PROJECT_NAME} PRIVATE
        OpenCL::OpenCL
        ${JSONCPP_LIBRARIES}
        #${CMAKE_THREAD_LIBS_INIT}
        Boost::filesystem
        ${VTK_LIBRARIES}
        stdc++
        m #needed to be added after implementing UpdatePrams() ??
)

set_target_properties(${PROJECT_NAME} PROPERTIES C_STANDARD 11
                                                 C_STANDARD_REQUIRED ON
                                                 C_EXTENSIONS OFF)

target_compile_definitions(${PROJECT_NAME} PRIVATE CL_TARGET_OPENCL_VERSION=300)
