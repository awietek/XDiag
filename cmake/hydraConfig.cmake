cmake_minimum_required(VERSION 3.19)

find_package(OpenMP)
  
set(hydra_known_comps static shared)
set(hydra_comp_static NO)
set(hydra_comp_shared NO)
foreach (hydra_comp IN LISTS ${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS)
    if (hydra_comp IN_LIST hydra_known_comps)
        set(hydra_comp_${hydra_comp} YES)
    else ()
        set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
            "hydra does not recognize component `${hydra_comp}`.")
        set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
        return()
    endif ()
endforeach ()

if (hydra_comp_static AND hydra_comp_shared)
    set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
        "hydra `static` and `shared` components are mutually exclusive.")
    set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
    return()
endif ()

set(hydra_static_targets "${CMAKE_CURRENT_LIST_DIR}/hydra-static-targets.cmake")
set(hydra_shared_targets "${CMAKE_CURRENT_LIST_DIR}/hydra-shared-targets.cmake")

macro(hydra_load_targets type)
    if (NOT EXISTS "${hydra_${type}_targets}")
        set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
            "hydra `${type}` libraries were requested but not found.")
        set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
        return()
    endif ()
    include("${hydra_${type}_targets}")
endmacro()

if (hydra_comp_static)
    hydra_load_targets(static)
elseif (hydra_comp_shared)
    hydra_load_targets(shared)
elseif (DEFINED hydra_SHARED_LIBS AND hydra_SHARED_LIBS)
    hydra_load_targets(shared)
elseif (DEFINED hydra_SHARED_LIBS AND NOT hydra_SHARED_LIBS)
    hydra_load_targets(static)
elseif (BUILD_SHARED_LIBS)
    if (EXISTS "${hydra_shared_targets}")
        hydra_load_targets(shared)
    else ()
        hydra_load_targets(static)
    endif ()
else ()
    if (EXISTS "${hydra_static_targets}")
        hydra_load_targets(static)
    else ()
        hydra_load_targets(shared)
    endif ()
endif ()