cmake_minimum_required(VERSION 3.20)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(LiftingLineTheory)

find_package(matplot++ CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

add_compile_definitions(${CMAKE_PROJECT_NAME} Airfoil_Folder="${CMAKE_SOURCE_DIR}/airfoil_data")
add_subdirectory(LLT)

add_executable(${CMAKE_PROJECT_NAME} main.cpp analysis_functions.cpp)


target_link_libraries(${CMAKE_PROJECT_NAME} Matplot++::matplot Eigen3::Eigen LLT)