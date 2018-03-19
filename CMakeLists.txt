goofit_add_directory()

#goofit_add_executable(D2K3_toy D2K3_toy.cu)
goofit_add_executable(D2KKK D2KKK.cpp)

# Add Meerkat
include(ExternalProject)
ExternalProject_Add(Meerkat
                    SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/meerkat"
                    URL "http://www.hepforge.org/archive/meerkat/meerkat-1.3.0.tar.gz"
                    CONFIGURE_COMMAND ""
                    BUILD_IN_SOURCE TRUE
                    INSTALL_COMMAND ""
                    )

add_library(MeerkatInt INTERFACE)
target_link_libraries(MeerkatInt INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/meerkat/lib/libMeerkat.so")
target_include_directories(MeerkatInt INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/meerkat/inc")
add_dependencies(MeerkatInt Meerkat)

target_link_libraries(D2KKK PUBLIC MeerkatInt)

# target_link_libraries(GooFit_Common INTERFACE ROOT::RooFit)
