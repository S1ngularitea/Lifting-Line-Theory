add_library(LLT LiftingLine.cpp)

find_package(Eigen3 CONFIG REQUIRED)

target_link_libraries(LLT Eigen3::Eigen)