cmake_minimum_required(VERSION 2.8)
project(SLIC6D)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -DBOOST_LOG_DYN_LINK")

find_package(OpenCV REQUIRED)
#message(STATUS "OpenCV_INCLUDE_DIRS = ${OpenCV_INCLUDE_DIRS}")
#message(STATUS "OpenCV_LIBS = ${OpenCV_LIBS}")
include_directories( ${OpenCV_INCLUDE_DIRS} )

find_package(Boost 1.40 COMPONENTS timer filesystem program_options log_setup log thread REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
add_definitions(-fext-numeric-literals)

add_subdirectory(slic6D_implementation)
set(ALL_MODULES slic_imp)

add_executable(run_slic run_slic.cpp)
target_link_libraries(run_slic ${ALL_MODULES} ${OpenCV_LIBS} ${Boost_LIBRARIES})

#add_executable(run_slic run_slic.cpp)
#target_link_libraries(run_slic ${OpenCV_LIBS})

add_executable(run_interpolation run_interpolation.cpp)
target_link_libraries(run_interpolation ${OpenCV_LIBS})
