#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/common.h>

#include <hydra/blocks/tj_symmetric/tj_symmetric_indexing.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms::tj_symmetric {

template <class bit_t, class coeff_t, class Filler>
void do_hopping_symmetric(BondList const &bonds, Couplings const &couplings,
                          indexing::tJSymmetricIndexing<bit_t> const &indexing,
                          Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;

  auto const &group_action = indexing.group_action();
  auto const &irrep = indexing.irrep();
  std::vector<coeff_t> bloch_factors;
  if constexpr (is_complex<coeff_t>()) {
    bloch_factors = irrep.characters();
  } else {
    bloch_factors = irrep.characters_real();
  }

  int n_sites = group_action.n_sites();

  auto clean_bonds =
      utils::clean_bondlist(bonds, couplings, {"HOP", "HOPUP", "HOPDN"}, 2);
  for (auto bond : clean_bonds) {

    std::string type = bond.type();
    std::string cpl = bond.coupling();

    assert(bond.size() == 2);
    int s1 = bond[0];
    int s2 = bond[1];
    if (s1 == s2) {
      lila::Log.err(
          "NotImplementedError: Hopping on two identical sites (tJ block)");
    }

    // Define hopping amplitude
    coeff_t t = 0;
    coeff_t tconj = 0;
    if constexpr (is_complex<coeff_t>()) {
      t = couplings[cpl];
      tconj = lila::conj(t);
    } else {
      t = lila::real(couplings[cpl]);
      tconj = lila::real(couplings[cpl]);
    }

    // Prepare bitmasks
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    int l = std::min(s1, s2);
    int u = std::max(s1, s2);
    bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);
    bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

    // Apply hoppings on dnspins
    if ((type == "HOP") || (type == "HOPDN")) {

      for (idx_t idx_up = 0; idx_up < indexing.n_reps_up(); ++idx_up) {
        bit_t ups = indexing.rep_up(idx_up);

        // skip if no hopping for dns works
        if (popcnt(ups & flipmask) != 0) {
          ++idx_up;
          continue;
        }

        idx_t up_offset = indexing.up_offset(idx_up);
        auto const &dnss = indexing.dns_for_up_rep(ups);
        auto [sym_begin, sym_end] = indexing.sym_limits_up(ups);

        // trivial stabilizer of ups -> dns have to be deposited
        if (sym_end - sym_begin == 1) {
          bit_t not_ups = (~ups) & sitesmask;
          idx_t idx_in = up_offset;
          for (bit_t dnsc : dnss) {
            bit_t dns = bitops::deposit(dnsc, not_ups);

            // If hopping is possible ...
            if (popcnt(dns & flipmask) == 1) {
              bit_t dns_flip = dns ^ flipmask;
              bit_t dns_flip_c = bitops::extract(dns_flip, not_ups);
              idx_t idx_out = up_offset + indexing.dnsc_index(dns_flip_c);

              // Complex conjugate t if necessary
              coeff_t val;
              if constexpr (is_complex<coeff_t>()) {
                val = -(gbit(dns, s1) ? t : tconj);
              } else {
                val = -t;
              } // Comment: norm is always 1.0 for trivial stabilizers

              // fill with correct fermi sign
              bool fermi = popcnt(dns & fermimask) & 1;
              fill(idx_out, idx_in, (fermi) ? -val : val);
            }

            ++idx_in;
          }
        }
        // non-trivial stabilizer of ups -> dns don't have to be deposited
        else {
          auto syms = indexing.syms_up(ups);
          auto const &norms = indexing.norms_for_up_rep(ups);

          idx_t idx_in = up_offset;
          idx_t idx_dns = 0;
          for (bit_t dns : dnss) {

            if (popcnt(dns & flipmask) == 1) {
              bit_t dns_flip = dns ^ flipmask;
              auto [idx_dns_flip, fermi_dn, sym] = indexing.index_dn_fermi_sym(
                  dns_flip, syms, dnss, fermimask);

              if (idx_dns_flip != invalid_index) {

                idx_t idx_out = up_offset + idx_dns_flip;
                bool fermi_up = indexing.fermi_bool_up(sym, ups);

                // Complex conjugate t if necessary
                coeff_t val;
                if constexpr (is_complex<coeff_t>()) {
                  val = -(gbit(dns, s1) ? t : tconj) * bloch_factors[sym] *
                        norms[idx_dns_flip] / norms[idx_dns];
                } else {
                  val = -t * bloch_factors[sym] * norms[idx_dns_flip] /
                        norms[idx_dns];
                }

                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
              }
            }

            ++idx_in;
            ++idx_dns;
          } // non-trivial stabilizer of ups

        } // for (auto [up, lower_upper]
      }
    } //     if ((hop.type() == "HOP") || (hop.type() == "HOPDN")) {

    // Apply hoppings on upspins
    if ((type == "HOP") || (type == "HOPUP")) {

      for (idx_t idx_up = 0; idx_up < indexing.n_reps_up(); ++idx_up) {
        bit_t ups = indexing.rep_up(idx_up);

        // continue if hopping not possible
        if (popcnt(ups & flipmask) != 1) {
          ++idx_up;
          continue;
        }

        bit_t ups_flip = ups ^ flipmask;
        idx_t idx_up_flip = indexing.index_up(ups_flip);
        bit_t ups_flip_rep = indexing.rep_up(idx_up_flip);
        bit_t not_ups_flip_rep = (~ups_flip_rep) & sitesmask;

        // Get limits, syms, and dns for ingoing ups
        idx_t up_offset_in = indexing.up_offset(idx_up);
        auto syms_up_in = indexing.syms_up(ups);
        auto const &dnss_in = indexing.dns_for_up_rep(ups);

        // Get limits, syms, and dns for outgoing ups
        idx_t up_offset_out = indexing.up_offset(idx_up_flip);
        auto syms_up_out = indexing.syms_up(ups_flip);
        auto const &dnss_out = indexing.dns_for_up_rep(ups_flip_rep);

        // auto const &dnss_out = indexing.dns_for_up_rep(ups_flip_rep);

        // Trivial stabilizer of target ups
        if (syms_up_out.size() == 1) {
          int sym = syms_up_out.front();

          // Complex conjugate t if necessary
          coeff_t prefac;
          if constexpr (is_complex<coeff_t>()) {
            prefac = -(gbit(ups, s1) ? t : tconj) * bloch_factors[sym];
          } else {
            prefac = -t * bloch_factors[sym];
          }

          // Fermi-sign of up spins
          bool fermi_up = (popcnt(ups & fermimask) & 1);
          fermi_up ^= indexing.fermi_bool_up(sym, ups_flip);

          // Origin ups trivial stabilizer -> dns need to be deposited
          if (syms_up_in.size() == 1) {
            idx_t idx_in = up_offset_in;
            idx_t idx_out = up_offset_out;
	    bit_t not_ups = (~ups) & sitesmask;
            for (bit_t dnsc : dnss_in) {
              bit_t dns = bitops::deposit(dnsc, not_ups);
              bool fermi_dn = popcnt(dns & fermimask) & 1;
              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -prefac : prefac);
              ++idx_in;
              ++idx_out;
            }
          }
          // Origin ups have stabilizer -> dns DONT need to be deposited
          else {
            auto const &norms_in = indexing.norms_for_up_rep(ups);
	    idx_t idx_dn = 0;
            idx_t idx_in = up_offset_in;
            for (bit_t dns : dnss_in) {
              auto [idx_dn_out, fermi_dn] = indexing.index_dn_fermi(
                  dns, sym, not_ups_flip_rep, fermimask);
              coeff_t val = prefac / norms_in[idx_dn];

              idx_t idx_out = up_offset_out + idx_dn_out;
              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
	      ++idx_dn;
              ++idx_in;
            }
          }
        }
        // Target ups have non-trivial stabilizer
        else {
          auto const &norms_out = indexing.norms_for_up_rep(ups_flip_rep);
          auto const &syms = syms_up_out;

          // Fix the bloch/prefactors
          std::vector<coeff_t> prefacs(irrep.size());
          if constexpr (is_complex<coeff_t>()) {
            for (int i = 0; i < (int)irrep.size(); ++i) {
              prefacs[i] = -(gbit(ups, s1) ? t : tconj) * irrep.character(i);
            }
          } else {
            for (int i = 0; i < (int)irrep.size(); ++i) {
              prefacs[i] = -t * lila::real(irrep.character(i));
            }
          }

          bool fermi_up_hop = (popcnt(ups & fermimask) & 1);

          // Origin ups trivial stabilizer -> dns need to be deposited
          if (syms_up_in.size() == 1) {
            idx_t idx_in = up_offset_in;
	    bit_t not_ups = (~ups) & sitesmask;
            for (bit_t dnsc : dnss_in) {
              bit_t dns = bitops::deposit(dnsc, not_ups);
              auto [idx_dn_out, fermi_dn, sym] =
                  indexing.index_dn_fermi_sym(dns, syms, dnss_out, fermimask);

              if (idx_dn_out != invalid_index) {
                idx_t idx_out = up_offset_out + idx_dn_out;
                bool fermi_up =
                    fermi_up_hop ^ indexing.fermi_bool_up(sym, ups_flip);
                coeff_t val = prefacs[sym] * norms_out[idx_dn_out];
                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
              }
              ++idx_in;
            }
          }
          // Origin ups non-trivial stabilizer -> dns DONT need to be deposited
          else {
            auto const &norms_in = indexing.norms_for_up_rep(ups);
            idx_t idx_in = up_offset_in;
	    idx_t idx_dn = 0;
            for (bit_t dns : dnss_in) {
              auto [idx_dn_out, fermi_dn, sym] =
                  indexing.index_dn_fermi_sym(dns, syms, dnss_out, fermimask);
              if (idx_dn_out != invalid_index) {
		idx_t idx_out = up_offset_out + idx_dn_out;
                bool fermi_up =
                    fermi_up_hop ^ indexing.fermi_bool_up(sym, ups_flip);
                coeff_t val =
                    prefacs[sym] * norms_out[idx_dn_out] / norms_in[idx_dn];
                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
              }
	      ++idx_dn;
              ++idx_in;
            }
          }
        } // Target ups have non-trivial stabilizer
      }   // loop over ups
    }     // type == "HOPUP"
  }       // for (auto bond : clean_bonds)
}
} // namespace hydra::terms::tj_symmetric
