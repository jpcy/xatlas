cmake_minimum_required(VERSION 3.10)
project(xatlas LANGUAGES CXX)

add_library(xatlas source/xatlas/xatlas.cpp)
add_library(xatlas::xatlas ALIAS xatlas)

target_include_directories(xatlas PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/source/xatlas>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS xatlas EXPORT xatlas-config 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(EXPORT xatlas-config NAMESPACE xatlas:: DESTINATION share/xatlas)
install(FILES source/xatlas/xatlas.h DESTINATION include)
