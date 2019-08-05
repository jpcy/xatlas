cmake_minimum_required(VERSION 3.0)

# Project is not versioned, so use date of last upstream commit.
project(xatlas LANGUAGES CXX VERSION "0.0.0")

set(XATLAS_SRC "xatlas.cpp")

set(XATLAS_HEADERS "xatlas.h")

add_library(xatlas ${XATLAS_SRC} ${XATLAS_HEADERS})

target_include_directories(
    xatlas
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>"
)

target_compile_definitions(xatlas PUBLIC "XA_DEBUG=$<CONFIG:Debug>")

# Introduce variables:
# * CMAKE_INSTALL_LIBDIR
# * CMAKE_INSTALL_BINDIR
# * CMAKE_INSTALL_INCLUDEDIR
include(GNUInstallDirs)

set(GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(VERSION_CONFIG "${GENERATED_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
set(PROJECT_CONFIG "${GENERATED_DIR}/${PROJECT_NAME}Config.cmake")

set(CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

set(TARGETS_EXPORT_NAME "${PROJECT_NAME}Targets")
set(NAMESPACE "${PROJECT_NAME}::")

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${VERSION_CONFIG}" COMPATIBILITY ExactVersion
)

configure_package_config_file(
    "cmake/Config.cmake.in"
    "${PROJECT_CONFIG}"
    INSTALL_DESTINATION "${CONFIG_INSTALL_DIR}"
)

install(
    TARGETS xatlas
    EXPORT "${TARGETS_EXPORT_NAME}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

install(
    FILES ${XATLAS_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)

install(
    FILES "${PROJECT_CONFIG}" "${VERSION_CONFIG}"
    DESTINATION "${CONFIG_INSTALL_DIR}"
)

install(
    EXPORT "${TARGETS_EXPORT_NAME}"
    NAMESPACE "${NAMESPACE}"
    DESTINATION "${CONFIG_INSTALL_DIR}"
)
