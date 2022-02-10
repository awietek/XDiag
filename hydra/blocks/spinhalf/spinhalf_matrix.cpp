#include "spinhalf_matrix.h"

#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/spinhalf/terms/spinhalf_terms.h>
#include <hydra/blocks/utils/block_utils.h>

#include <hydra/operators/operator_qns.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
lila::Matrix<coeff_t>
MatrixGen(BondList const &bonds, Couplings const &couplings,
          Spinhalf<bit_t> const &block_in, Spinhalf<bit_t> const &block_out) {

  int n_up_out = utils::spinhalf_nup(bonds, couplings, block_in);
  if (n_up_out != block_out.n_up())
    lila::Log.err("Incompatible n_up in MatrixGen: {} != {}", n_up_out,
                  block_out.n_up());

  utils::check_operator_works_with<coeff_t>(bonds, couplings,
                                            "spinhalf_matrix");
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<coeff_t>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  auto const &indexing_out = block_out.indexing();

  if (block_in == block_out) {
    terms::spinhalf_ising<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
    terms::spinhalf_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                             fill);
    terms::spinhalf_sz<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  }
  terms::spinhalf_spsm<bit_t, coeff_t>(bonds, couplings, indexing_in,
				       indexing_out, fill, "S+");
  terms::spinhalf_spsm<bit_t, coeff_t>(bonds, couplings, indexing_in,
				       indexing_out, fill, "S-");

  return mat;
}

template lila::Matrix<double>
MatrixGen<uint16_t, double>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint16_t> const &block_in,
                            Spinhalf<uint16_t> const &block_out);
template lila::Matrix<double>
MatrixGen<uint32_t, double>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint32_t> const &block_in,
                            Spinhalf<uint32_t> const &block_out);
template lila::Matrix<double>
MatrixGen<uint64_t, double>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint64_t> const &block_in,
                            Spinhalf<uint64_t> const &block_out);

template lila::Matrix<complex>
MatrixGen<uint16_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Spinhalf<uint16_t> const &block_in,
                             Spinhalf<uint16_t> const &block_out);
template lila::Matrix<complex>
MatrixGen<uint32_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Spinhalf<uint32_t> const &block_in,
                             Spinhalf<uint32_t> const &block_out);
template lila::Matrix<complex>
MatrixGen<uint64_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Spinhalf<uint64_t> const &block_in,
                             Spinhalf<uint64_t> const &block_out);

} // namespace hydra
