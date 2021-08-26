#include "electron_symmetric_matrix.h"

#include <hydra/utils/bitops.h>

#include <hydra/models/electron_symmetric/terms/electron_symmetric_exchange.h>
#include <hydra/models/electron_symmetric/terms/electron_symmetric_hopping.h>
#include <hydra/models/electron_symmetric/terms/electron_symmetric_ising.h>
#include <hydra/models/electron_symmetric/terms/electron_symmetric_u.h>

namespace hydra {

template <class bit_t, class GroupAction>
lila::Matrix<double>
MatrixReal(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, GroupAction> const &block_in,
           ElectronSymmetric<bit_t, GroupAction> const &block_out) {
  assert(block_in == block_out); // only temporary

  if (is_complex(bonds))
    lila::Log.err("Cannot construct real matrix from complex bonds!");
  if (is_complex(couplings))
    lila::Log.err("Cannot construct real matrix from complex couplings!");
  if (is_complex(block_in.irrep()) || is_complex(block_out.irrep()))
    lila::Log.err("Cannot construct real matrix from complex representation!");

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<double>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, double val) {
    mat(idx_out, idx_in) += val;
  };

  electron::do_U_symmetric(couplings, block_in, fill);
  electron::do_hopping_symmetric<bit_t, double>(bonds, couplings, block_in,
                                                fill);
  electron::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  electron::do_exchange_symmetric<bit_t, double>(bonds, couplings, block_in,
                                                 fill);
  return mat;
}

template lila::Matrix<double> MatrixReal<uint16, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16, PermutationGroupAction> const &block_in,
    ElectronSymmetric<uint16, PermutationGroupAction> const &block_out);
template lila::Matrix<double> MatrixReal<uint32, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32, PermutationGroupAction> const &block_in,
    ElectronSymmetric<uint32, PermutationGroupAction> const &block_out);
template lila::Matrix<double> MatrixReal<uint64, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64, PermutationGroupAction> const &block_in,
    ElectronSymmetric<uint64, PermutationGroupAction> const &block_out);

template <class bit_t, class GroupAction>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, GroupAction> const &block_in,
           ElectronSymmetric<bit_t, GroupAction> const &block_out) {
  assert(block_in == block_out); // only temporary

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<complex>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, complex val) {
    mat(idx_out, idx_in) += val;
  };

  electron::do_U_symmetric(couplings, block_in, fill);
  electron::do_hopping_symmetric<bit_t, complex>(bonds, couplings, block_in,
                                                 fill);
  electron::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  electron::do_exchange_symmetric<bit_t, complex>(bonds, couplings, block_in,
                                                  fill);
  return mat;
}

template lila::Matrix<complex> MatrixCplx<uint16, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16, PermutationGroupAction> const &block_in,
    ElectronSymmetric<uint16, PermutationGroupAction> const &block_out);
template lila::Matrix<complex> MatrixCplx<uint32, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32, PermutationGroupAction> const &block_in,
    ElectronSymmetric<uint32, PermutationGroupAction> const &block_out);
template lila::Matrix<complex> MatrixCplx<uint64, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64, PermutationGroupAction> const &block_in,
    ElectronSymmetric<uint64, PermutationGroupAction> const &block_out);

} // namespace hydra