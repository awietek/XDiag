cmake_minimum_required(VERSION 3.19)

project(
  hb_chain_dynamical_sf_fulled
)

add_executable(main main.cpp)

# set xdiag compile options
find_package(xdiag REQUIRED HINTS "../../../install")
target_link_libraries(main PRIVATE xdiag::xdiag)
