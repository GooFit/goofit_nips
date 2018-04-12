cmake_minimum_required(VERSION 3.8)

goofit_add_package(NIPS)

goofit_add_executable(D2KKK D2KKK.cu)



if(GOOFIT_DEVICE STREQUAL CUDA)
    target_compile_options(D2KKK PUBLIC --expt-extended-lambda)
endif()
