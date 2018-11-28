cmake_minimum_required(VERSION 3.1.0)

enable_language(Fortran)
project(gauss-glmnet)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set_source_files_properties(glmnet.f PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -ffixed-form")

add_library(glmnet SHARED glmnet.cpp glmnet5dp.f)

