cmake_minimum_required(VERSION 2.8)
project(SLIC6D)

find_package(OpenCV REQUIRED)
message(STATUS "OpenCV_INCLUDE_DIRS = ${OpenCV_INCLUDE_DIRS}")
message(STATUS "OpenCV_LIBS = ${OpenCV_LIBS}")
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_subdirectory(slic6D_implementation)
set(ALL_MODULES slic_imp)

add_executable(run_slic run_slic.cpp)
target_link_libraries(run_slic ${ALL_MODULES} ${OpenCV_LIBS})

#add_executable(run_slic run_slic.cpp)
#target_link_libraries(run_slic ${OpenCV_LIBS})

add_executable(run_interpolation run_interpolation.cpp)
target_link_libraries(run_interpolation ${OpenCV_LIBS})
