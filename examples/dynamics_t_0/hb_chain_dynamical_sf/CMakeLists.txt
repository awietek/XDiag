cmake_minimum_required(VERSION 3.19)

project(
  hb_chain_dynamical_sf
)

add_executable(main main.cpp)
find_package(xdiag REQUIRED HINTS "../../../install")
target_link_libraries(main PUBLIC xdiag::xdiag)
