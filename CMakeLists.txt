cmake_minimum_required(VERSION 3.8)

goofit_add_directory()

goofit_add_executable(D2KKK D2KKK.cu)

## Add Meerkat
if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/meerkat/lib/libMeerkat.so")
    include(ExternalProject)
    ExternalProject_Add(Meerkat
                        SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/meerkat"
                        URL "http://www.hepforge.org/archive/meerkat/meerkat-1.3.0.tar.gz"
                        CONFIGURE_COMMAND ""
                        BUILD_IN_SOURCE TRUE
                        INSTALL_COMMAND ""
                        )
endif()

add_library(MeerkatInt INTERFACE)
target_link_libraries(MeerkatInt INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/meerkat/lib/libMeerkat.so")
target_link_libraries(MeerkatInt INTERFACE ROOT::RooFit ROOT::RooFitCore)
target_include_directories(MeerkatInt INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/meerkat/inc")

if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/meerkat/lib/libMeerkat.so")
    add_dependencies(MeerkatInt Meerkat)
endif()

target_link_libraries(D2KKK PUBLIC MeerkatInt)

if(GOOFIT_DEVICE STREQUAL CUDA)
    target_compile_options(D2KKK PUBLIC --expt-extended-lambda)
endif()
