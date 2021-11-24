#pragma once

#include <hydra/common.h>

#include <hydra/blocks/utils/block_utils.h>

#include <hydra/bitops/bitops.h>

#include <hydra/indexing/tj/tj_symmetric_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms::tj_symmetric {

template <typename bit_t, typename coeff_t, class Filler>
void do_exchange_symmetric(BondList const &bonds, Couplings const &couplings,
                           indexing::tJSymmetricIndexing<bit_t> const &indexing,
                           Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;

  auto const &group_action = indexing.group_action();
  auto const &irrep = indexing.irrep();
  int n_sites = group_action.n_sites();

  auto clean_bonds = utils::clean_bondlist(
      bonds, couplings,
      {"HEISENBERG", "HB", "EXCHANGE", "TJHEISENBERG", "TJHB"}, 2);
  for (auto bond : clean_bonds) {

    std::string cpl = bond.coupling();

    auto [J, J_conj] = utils::get_coupling_and_conj<coeff_t>(couplings, cpl);

    utils::check_sites_disjoint(bond);
    int s1 = bond[0];
    int s2 = bond[1];

    // Prepare bitmasks
    bit_t s1mask = ((bit_t)1 << s1);
    bit_t s2mask = ((bit_t)1 << s2);
    bit_t flipmask = s1mask | s2mask;
    int l = std::min(s1, s2);
    int u = std::max(s1, s2);
    bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);
    bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

    // // DEBUG PRINT states
    // lila::Log("\n\n\n");
    // idx_t idx = 0;
    // for (idx_t idx_up = 0; idx_up < indexing.n_rep_ups(); ++idx_up) {
    //   bit_t ups = indexing.rep_up(idx_up);
    //   auto syms = indexing.syms_up(ups);
    //   auto const &dnss = indexing.dns_for_up_rep(ups);
    //   if (syms.size() == 1) {
    //     bit_t not_ups = (~ups) & sitesmask;

    //     for (bit_t dnsc : dnss) {
    //       bit_t dns = bitops::deposit(dnsc, not_ups);
    //       lila::Log("{}: {};{}", idx, BSTR(ups), BSTR(dns));
    //       ++idx;
    //     }
    //   } else {
    //     for (bit_t dns : dnss) {
    //       lila::Log("{}: {};{}", idx, BSTR(ups), BSTR(dns));
    //       ++idx;
    //     }
    //   }
    // }

    // Loop over all up configurations
    for (idx_t idx_ups = 0; idx_ups < indexing.n_rep_ups(); ++idx_ups) {
      bit_t ups = indexing.rep_ups(idx_ups);
      bit_t not_ups = (~ups) & sitesmask;

      // Exchange gives zero continue
      if (gbit(ups, s1) == gbit(ups, s2)) {
        continue;
      }

      // Set the correct prefactor
      coeff_t Jhalf;
      if constexpr (is_complex<coeff_t>()) {
        Jhalf = (gbit(ups, s1)) ? J / 2. : J_conj / 2.;
      } else {
        Jhalf = J / 2.;
      }

      // Get the flip masks for up and down spins
      bit_t dns_mask = ((ups & flipmask) == s1mask) ? s2mask : s1mask;

      // Compute index and rep of flipped ups
      bit_t ups_flip = ups ^ flipmask;
      idx_t idx_ups_flip = indexing.index_ups(ups_flip);
      bit_t ups_flip_rep = indexing.rep_ups(idx_ups_flip);
      bit_t not_ups_flip_rep = (~ups_flip_rep) & sitesmask;

      // Get limits, syms, and dns for ingoing ups
      idx_t up_offset_in = indexing.ups_offset(idx_ups);
      auto syms_up_in = indexing.syms_ups(ups);
      auto dnss_in = indexing.dns_for_ups_rep(ups);

      // Get limits, syms, and dns for outgoing ups
      idx_t up_offset_out = indexing.ups_offset(idx_ups_flip);
      auto syms_up_out = indexing.syms_ups(ups_flip);
      auto dnss_out = indexing.dns_for_ups_rep(ups_flip_rep);

      // Target ups have trivial stabilizer
      if (syms_up_out.size() == 1) {
        int sym = syms_up_out[0];

        // Fix the bloch factor
        coeff_t prefac;
        if constexpr (is_complex<coeff_t>()) {
          prefac = -Jhalf * irrep.character(sym);
        } else {
          prefac = -Jhalf * lila::real(irrep.character(sym));
        }

        // Fermi-sign of up spins
        bool fermi_up = (popcnt(ups & fermimask) & 1);
        fermi_up ^= indexing.fermi_bool_ups(sym, ups_flip);

        // Origin ups trivial stabilizer -> dns need to be deposited
        if (syms_up_in.size() == 1) {
          idx_t idx_dn = 0;
          for (bit_t dnsc : dnss_in) {
            bit_t dns = bitops::deposit(dnsc, not_ups);

            if ((dns & flipmask) == dns_mask) { // If  dns can be raised
              bit_t dns_flip = dns ^ flipmask;

              // lila::Log("-----------------------------------------");
              // lila::Log("CASE: origin-stab-FALSE, target-stab-FALSE");
              // lila::Log("origin-sym: {}", sym);
              // lila::Log("from: {};{}", BSTR(ups), BSTR(dns));
              // lila::Log("mask: {} {}", BSTR(flipmask), BSTR(flipmask));
              // lila::Log("to  : {};{}", BSTR(ups_flip), BSTR(dns_flip));
              // bit_t dns_flip_rep = group_action.apply(sym, dns_flip);
              // lila::Log("rep : {};{}", BSTR(ups_flip_rep),
              // BSTR(dns_flip_rep));

              auto [idx_dn_flip, fermi_dn] = indexing.index_dns_fermi(
                  dns_flip, sym, not_ups_flip_rep, fermimask);

              // lila::Log("fermi_up: {}, fermi_dn: {}", fermi_up, fermi_dn);

              coeff_t val = prefac; // norms are both 1.0 in this case
              idx_t idx_in = up_offset_in + idx_dn;
              idx_t idx_out = up_offset_out + idx_dn_flip;
              // lila::Log("idx_dn: {}, idx_dn_flip: {}", idx_dn, idx_dn_flip);

              // lila::Log("off_in: {}, off_out: {}", up_offset_in,
              // up_offset_out); lila::Log("idx_in: {}, idx_out: {}", idx_in,
              // idx_out); lila::Log("val: {}", lila::real(val));
              // lila::Log("fill: {}", lila::real((fermi_up ^ fermi_dn) ? -val :
              // val));

              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
            }

            ++idx_dn;
          }
        }
        // Origin ups non-trivial stabilizer -> dns DONT need to be deposited
        else {

          auto norms_in = indexing.norms_for_ups_rep(ups);

          idx_t idx_dn = 0;
          for (bit_t dns : dnss_in) {

            if ((dns & flipmask) == dns_mask) { // If  dns can be raised
              bit_t dns_flip = dns ^ flipmask;

              // lila::Log("-----------------------------------------");
              // lila::Log("CASE: origin-stab_TRUE, target-stab-FALSE");
              // lila::Log("origin-stab:");
              // for (auto sym : syms_up_in) {
              //   std::cout << sym << " ";
              // }
              // lila::Log("");

              // lila::Log("from: {};{}", BSTR(ups), BSTR(dns));
              // lila::Log("mask: {} {}", BSTR(flipmask), BSTR(flipmask));
              // lila::Log("to  : {},{}", BSTR(ups_flip), BSTR(dns_flip));
              // bit_t dns_flip_rep = group_action.apply(sym, dns_flip);
              // lila::Log("rep : {},{}", BSTR(ups_flip_rep),
              // BSTR(dns_flip_rep));

              auto [idx_dn_flip, fermi_dn] = indexing.index_dns_fermi(
                  dns_flip, sym, not_ups_flip_rep, fermimask);

              // lila::Log("idx_dn: {}, idx_dn_flip: {}", idx_dn, idx_dn_flip);
              // lila::Log("fermi_up: {}, fermi_dn: {}", fermi_up, fermi_dn);

              coeff_t val = prefac / norms_in[idx_dn];

              idx_t idx_in = up_offset_in + idx_dn;
              idx_t idx_out = up_offset_out + idx_dn_flip;
              // lila::Log("idx_in: {}, idx_out: {}", idx_in, idx_out);
              // lila::Log("val: {}", lila::real(val));
              // lila::Log("fill: {}", lila::real((fermi_up ^ fermi_dn) ? -val :
              // val));
              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
            }
            ++idx_dn;
          }
        }
      }
      // Target ups have non-trivial stabilizer
      else {
        auto norms_out = indexing.norms_for_ups_rep(ups_flip_rep);
        auto syms = syms_up_out;

        // Fix the bloch/prefactors
        std::vector<coeff_t> prefacs(irrep.size());
        if constexpr (is_complex<coeff_t>()) {
          for (int i = 0; i < (int)irrep.size(); ++i) {
            prefacs[i] = -irrep.character(i) * Jhalf;
          }
        } else {
          for (int i = 0; i < (int)irrep.size(); ++i) {
            prefacs[i] = -lila::real(irrep.character(i)) * Jhalf;
          }
        }

        bool fermi_up_hop = (popcnt(ups & fermimask) & 1);

        // Origin ups trivial stabilizer -> dns need to be deposited
        if (syms_up_in.size() == 1) {
          idx_t idx_dn = 0;
          for (bit_t dnsc : dnss_in) {
            bit_t dns = bitops::deposit(dnsc, not_ups);

            if ((dns & flipmask) == dns_mask) { // If  dns can be raised
              bit_t dns_flip = dns ^ flipmask;

              // lila::Log("-----------------------------------------");
              // lila::Log("CASE: origin-stab-FALSE, target-stab-TRUE");
              // lila::Log("target-stab:");
              // for (auto sym : syms_up_out) {
              //   std::cout << sym << " ";
              // }
              // lila::Log("");

              // lila::Log("from: {};{}", BSTR(ups), BSTR(dns));
              // lila::Log("mask: {} {}", BSTR(flipmask), BSTR(flipmask));
              // lila::Log("to  : {},{}", BSTR(ups_flip), BSTR(dns_flip));
              // bit_t dns_flip_rep = symmetries::representative_subset(
              //     dns_flip, group_action, syms);
              // lila::Log("rep : {};{}", BSTR(ups_flip_rep),
              // BSTR(dns_flip_rep));
              auto [idx_dn_flip, fermi_dn, sym] = indexing.index_dns_fermi_sym(
                  dns_flip, syms, dnss_out, fermimask);

              // lila::Log("idx_dn: {}, idx_dn_flip: {}", idx_dn, idx_dn_flip);

              if (idx_dn_flip != invalid_index) {
                bool fermi_up =
                    fermi_up_hop ^ indexing.fermi_bool_ups(sym, ups_flip);

                // lila::Log("fermi_up: {}, fermi_dn: {}", fermi_up, fermi_dn);

                coeff_t val = prefacs[sym] * norms_out[idx_dn_flip];
                idx_t idx_in = up_offset_in + idx_dn;
                idx_t idx_out = up_offset_out + idx_dn_flip;
                // lila::Log("idx_in: {}, idx_out: {}", idx_in, idx_out);
                // lila::Log("val: {}", lila::real(val));

                // lila::Log("fill: {}", lila::real((fermi_up ^ fermi_dn) ? -val
                // : val));

                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
              }
            }

            ++idx_dn;
          }
        }
        // Origin ups non-trivial stabilizer -> dns DONT need to be deposited
        else {
          auto norms_in = indexing.norms_for_ups_rep(ups);

          idx_t idx_dn = 0;
          for (bit_t dns : dnss_in) {
            if ((dns & flipmask) == dns_mask) { // If  dns can be raised
              bit_t dns_flip = dns ^ flipmask;

              // lila::Log("-----------------------------------------");
              // lila::Log("CASE: origin-stab-TRUE, target-stab-TRUE");
              // lila::Log("origin-stab:");
              // for (auto sym : syms_up_in) {
              //   std::cout << sym << " ";
              // }
              // lila::Log("");
              // lila::Log("target-stab:");
              // for (auto sym : syms_up_out) {
              //   std::cout << sym << " ";
              // }
              // lila::Log("");

              // lila::Log("from: {};{}", BSTR(ups), BSTR(dns));
              // lila::Log("mask: {} {}", BSTR(flipmask), BSTR(flipmask));
              // lila::Log("to  : {};{}", BSTR(ups_flip), BSTR(dns_flip));
              // bit_t dns_flip_rep = symmetries::representative_subset(
              //     dns_flip, group_action, syms);
              // lila::Log("rep : {};{}", BSTR(ups_flip_rep),
              // BSTR(dns_flip_rep));
              auto [idx_dn_flip, fermi_dn, sym] = indexing.index_dns_fermi_sym(
                  dns_flip, syms, dnss_out, fermimask);

              // lila::Log("idx_dn: {}, idx_dn_flip: {}", idx_dn, idx_dn_flip);

              if (idx_dn_flip != invalid_index) {

                bool fermi_up =
                    fermi_up_hop ^ indexing.fermi_bool_ups(sym, ups_flip);

                // lila::Log("fermi_up: {}, fermi_dn: {}", fermi_up, fermi_dn);

                coeff_t val =
                    prefacs[sym] * norms_out[idx_dn_flip] / norms_in[idx_dn];
                idx_t idx_in = up_offset_in + idx_dn;
                idx_t idx_out = up_offset_out + idx_dn_flip;

                // lila::Log("idx_in: {}, idx_out: {}", idx_in, idx_out);
                // lila::Log("val: {}", lila::real(val));
                // lila::Log("fill: {}", lila::real((fermi_up ^ fermi_dn) ? -val
                // : val));
                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
              }
            }

            ++idx_dn;
          }
        }
        // for (idx_t idx_up = 0; ...
      } // for (auto bond : clean_bonds)
    }
  }
}
} // namespace hydra::terms::tj_symmetric
