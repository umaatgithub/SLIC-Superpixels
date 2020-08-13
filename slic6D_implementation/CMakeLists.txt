find_package(OpenCV REQUIRED)
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_library(slic_imp slic6d.cpp slic6d.h)

target_link_libraries(slic_imp ${OpenCV_LIBS})

