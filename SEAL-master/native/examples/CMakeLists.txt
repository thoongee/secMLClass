# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.10)

#project(SEALExamples VERSION 3.2.0 LANGUAGES CXX)
project(SEALMLClass VERSION 3.2.0 LANGUAGES CXX)

# Executable will be in ../bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin)

# set C++ version (jina add)
set(CMAKE_CXX_STANDARD 11)

# add because of std::experimental::filesystem::v1::__cxx11 error
link_libraries(stdc++fs)

#add_executable(sealexamples examples.cpp)
add_executable(secMLClass secMLClass.cpp strutils.cpp randgen.cpp NaiveBayesBreastCancer.cpp MultinomialNB_Email.cpp)

# Import Microsoft SEAL
find_package(SEAL 3.2.0 EXACT REQUIRED)

# Link Microsoft SEAL
#target_link_libraries(sealexamples SEAL::seal)
target_link_libraries(secMLClass SEAL::seal) 
