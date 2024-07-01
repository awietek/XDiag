#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/bondlist.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

struct eigvals_lanczos_result_t {
  arma::vec alphas;
  arma::vec betas;
  arma::vec eigenvalues;
  int64_t niterations;
  std::string criterion;
};

eigvals_lanczos_result_t
eigvals_lanczos(BondList const &bonds, block_variant_t const &block,
                int64_t neigvals = 1, double precision = 1e-12,
                int64_t max_iterations = 1000, bool force_complex = false,
                double deflation_tol = 1e-7, int64_t random_seed = 42);

eigvals_lanczos_result_t
eigvals_lanczos(BondList const &bonds, block_variant_t const &block,
                State &state0,
                int64_t neigvals = 1, double precision = 1e-12,
                int64_t max_iterations = 1000, bool force_complex = false,
                double deflation_tol = 1e-7);

} // namespace xdiag
