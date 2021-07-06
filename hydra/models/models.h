#pragma once

#include <hydra/common.h>

namespace hydra {

template <class bit_t> class Spinhalf;
template <class bit_t> class SpinhalfMPI;

template <class bit_t> class tJ;

template <class bit_t> class Electron;
template <class bit_t, class SymmetryGroup> class ElectronSymmetric;

template <class Block> struct is_mpi_block_t {
  static constexpr bool value = false;
};
template <> struct is_mpi_block_t<SpinhalfMPI<uint16>> {
  static constexpr bool value = true;
};
template <> struct is_mpi_block_t<SpinhalfMPI<uint32>> {
  static constexpr bool value = true;
};
template <> struct is_mpi_block_t<SpinhalfMPI<uint64>> {
  static constexpr bool value = true;
};

template <class Block>
constexpr bool is_mpi_block = is_mpi_block_t<Block>::value;

} // namespace hydra
