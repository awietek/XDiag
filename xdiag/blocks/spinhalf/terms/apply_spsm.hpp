#pragma once

#include <functional>
#include <string>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/blocks/spinhalf/terms/apply_term_offdiag_no_sym.hpp>
#include <xdiag/blocks/spinhalf/terms/apply_term_offdiag_sym.hpp>
#include <xdiag/common.hpp>

namespace xdiag::spinhalf {

// S+ or S- term: J S^+_i   OR   J S^-_i

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_spsm(Bond const &bond, BasisIn &&basis_in, BasisOut &&basis_out,
                Fill &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert((bond.type() == "S+") || (bond.type() == "S-"));
  assert(bond.size() == 1);

  coeff_t J = bond.coupling<coeff_t>();
  int64_t s = bond[0];
  bit_t mask = ((bit_t)1 << s);

  // Define actions of bond
  std::function<bool(bit_t)> non_zero_term;
  std::function<std::pair<bit_t, coeff_t>(bit_t)> term_action;

  if (bond.type() == "S+") {
    non_zero_term = [&mask](bit_t spins) { return !(spins & mask); };
    term_action = [&mask, &J](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bit_t spins_flip = spins | mask;
      return {spins_flip, J};
    };
  } else { // bond.type() == "S-"
    non_zero_term = [&mask](bit_t spins) { return spins & mask; };
    term_action = [&mask, &J](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bit_t spins_flip = spins ^ mask;
      return {spins_flip, J};
    };
  }

  // Dispatch either symmetric of unsymmetric term application
  if constexpr (symmetric) {
    spinhalf::apply_term_offdiag_sym<bit_t, coeff_t>(
        basis_in, basis_out, non_zero_term, term_action, fill);
  } else {
    spinhalf::apply_term_offdiag_no_sym<bit_t, coeff_t>(
        basis_in, basis_out, non_zero_term, term_action, fill);
  }
}

} // namespace xdiag::spinhalf
