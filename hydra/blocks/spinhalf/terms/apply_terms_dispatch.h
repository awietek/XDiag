#pragma once

#include <functional>

#include <hydra/basis/basis.h>
#include <hydra/blocks/spinhalf/terms/apply_terms.h>
#include <hydra/common.h>
#include <hydra/utils/logger.h>

namespace hydra {

template <typename coeff_t, class fill_f>
void apply_terms_dispatch(BondList const &bonds, Spinhalf const &block_in,
                          Spinhalf const &block_out, fill_f fill) try {
  using namespace spinhalf;
  using namespace basis::spinhalf;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      overload{
          // uint16_t
          [&](BasisSz<uint16_t> const &idx_in,
              BasisSz<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisNoSz<uint16_t> const &idx_in,
              BasisNoSz<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricSz<uint16_t> const &idx_in,
              BasisSymmetricSz<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNoSz<uint16_t> const &idx_in,
              BasisSymmetricNoSz<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint16_t, 1> const &idx_in,
              BasisSublattice<uint16_t, 1> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint16_t, 2> const &idx_in,
              BasisSublattice<uint16_t, 2> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint16_t, 3> const &idx_in,
              BasisSublattice<uint16_t, 3> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint16_t, 4> const &idx_in,
              BasisSublattice<uint16_t, 4> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint16_t, 5> const &idx_in,
              BasisSublattice<uint16_t, 5> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },

          // uint32_t
          [&](BasisSz<uint32_t> const &idx_in,
              BasisSz<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisNoSz<uint32_t> const &idx_in,
              BasisNoSz<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricSz<uint32_t> const &idx_in,
              BasisSymmetricSz<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNoSz<uint32_t> const &idx_in,
              BasisSymmetricNoSz<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint32_t, 1> const &idx_in,
              BasisSublattice<uint32_t, 1> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint32_t, 2> const &idx_in,
              BasisSublattice<uint32_t, 2> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint32_t, 3> const &idx_in,
              BasisSublattice<uint32_t, 3> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint32_t, 4> const &idx_in,
              BasisSublattice<uint32_t, 4> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint32_t, 5> const &idx_in,
              BasisSublattice<uint32_t, 5> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },

          // uint64_t
          [&](BasisSz<uint64_t> const &idx_in,
              BasisSz<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisNoSz<uint64_t> const &idx_in,
              BasisNoSz<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricSz<uint64_t> const &idx_in,
              BasisSymmetricSz<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNoSz<uint64_t> const &idx_in,
              BasisSymmetricNoSz<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint64_t, 1> const &idx_in,
              BasisSublattice<uint64_t, 1> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint64_t, 2> const &idx_in,
              BasisSublattice<uint64_t, 2> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint64_t, 3> const &idx_in,
              BasisSublattice<uint64_t, 3> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint64_t, 4> const &idx_in,
              BasisSublattice<uint64_t, 4> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSublattice<uint64_t, 5> const &idx_in,
              BasisSublattice<uint64_t, 5> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },

          [&](auto const &idx_in, auto const &idx_out) {
            HydraThrow(std::logic_error, "Invalid basis or combination of bases");
            (void)idx_in;
            (void)idx_out;
          }},
      basis_in, basis_out);
} catch (...) {
  HydraRethrow("Unable to apply terms on Spinhalf block");
}

} // namespace hydra
