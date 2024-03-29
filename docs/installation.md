---
title: Installation
---

## C++ Compilation

Using hydra with C++ is a two-step process. First the `hydra` library needs
to be compiled and installed. Therafter, application codes are compiled
in a second step. Here we explain how to compile the library.

### Prerequisites

* A C++ compiler that supports C++17 (`g++`, `clang`, or Intel's `icpx`)
* [git](https://git-scm.com/) version control system
* [CMake](https://cmake.org/) build system generator 
* A linear algebra backend (Blas/Lapack, IntelMKL or Accelerate on OSX)
* **optional** [HDF5](https://www.hdfgroup.org/solutions/hdf5/), [OpenMP](https://www.openmp.org/)
* **optional** [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) for the distributed library

### Basic Compilation

- **Download** the source code using [git](https://git-scm.com/)
  ```bash
  cd /path/where/hydra/should/be
  git clone https://github.com/awietek/hydra.git
  ```

- **Compile the default library**
  ``` bash
  cd hydra
  cmake -S . -B build
  cmake --build build
  cmake --install build
  ```
  By default, the library is now installed in the subdirectory `install`.

- **Compile the distributed library**

    To use the distributed computing features of `hydra`, the distributed
    library has to be built which requires [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface).
    ``` bash
    cd hydra
    cmake -S . -B build -D HYDRA_DISTRIBUTED=On
    cmake --build build
    cmake --install build
    ```

### Advanced Compilation

- **Parallel compilation**
    To speed up the compilation process, the build step can be performed in parallel using the `-j` flag

    ```bash
    cmake --build build -j
    ```

- **Listing compile options**

    The available compilation options can be displayed using
    ``` bash
    cmake -L .
    ```

- **Choosing a certain compiler**

    The compiler (e.g. `icpx`) can be specified using
    ``` bash
    cmake -S . -B build -D CMAKE_CXX_COMPILER=icpx
    ```

    !!! warning 

        If the `hydra` library is compiled with a certain compiler, it is
        advisable to also compile the application codes with the same compiler.

- **Setting the install path**

    In the installation step, the install directory can be set in the following way
    ```bash
    cmake --install build --prefix /my/install/prefix
    ```

- **Disabling HDF5/OpenMP**

    To disable support for [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
    or [OpenMP](https://www.openmp.org/) support, use
    ```bash
    cmake -S . -B build -D HYDRA_DISABLE_OPENMP=On -D HYDRA_DISABLE_HDF5=On
    ```
    
- **Building and running tests**

    To compile and run the testing programs, use
    ``` bash
    cmake -S . -B build -D BUILD_TESTING=On
    cmake --build build
    build/tests/tests
    ```