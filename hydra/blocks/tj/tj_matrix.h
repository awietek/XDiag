#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/common.h>

#include <hydra/blocks/blocks.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
arma::Mat<coeff_t> MatrixGen(BondList const &bonds, Couplings const &couplings,
                             tJ<bit_t> const &block_in,
                             tJ<bit_t> const &block_out);

template <typename bit_t>
arma::Mat<double> MatrixReal(BondList const &bonds, Couplings const &couplings,
                             tJ<bit_t> const &block_in,
                             tJ<bit_t> const &block_out) {
  return MatrixGen<bit_t, double>(bonds, couplings, block_in, block_out);
}
template <typename bit_t>
arma::Mat<complex> MatrixCplx(BondList const &bonds, Couplings const &couplings,
                              tJ<bit_t> const &block_in,
                              tJ<bit_t> const &block_out) {
  return MatrixGen<bit_t, complex>(bonds, couplings, block_in, block_out);
}

} // namespace hydra
