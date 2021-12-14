#pragma once

#include <map>
#include <vector>

#include <lila/external/gsl/span>

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/permutation_group_lookup.h>
#include <hydra/symmetries/representation.h>
#include <hydra/symmetries/symmetry_operations.h>
#include <hydra/symmetries/representative_list.h>
#include <hydra/symmetries/fermi_bool_table.h>

namespace hydra::indexing {

template <typename bit_t> class tJSymmetricIndexing {
public:
  using span_size_t = gsl::span<int const>::size_type;

  tJSymmetricIndexing(int n_sites, int nup, int ndn,
                      PermutationGroup permutation_group, Representation irrep);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }
  PermutationGroupLookup<bit_t> const &group_action() const {
    return group_action_;
  }
  Representation const &irrep() const { return irrep_; }
  idx_t size() const { return size_; }

  idx_t n_rep_ups() const { return reps_up_.size(); }
  bit_t rep_ups(idx_t idx_up) const { return reps_up_[idx_up]; }
  idx_t ups_offset(idx_t idx_up) const { return ups_offset_[idx_up]; }

  // index and fermi sign for dns with trivial stabilizer
  std::pair<idx_t, bool> index_dns_fermi(bit_t dns, int sym,
                                         bit_t not_ups) const {
    bit_t dns_rep = group_action_.apply(sym, dns);
    bit_t dns_rep_c = bitops::extract(dns_rep, not_ups);
    idx_t idx_dns_rep = lintable_dnsc_.index(dns_rep_c);
    bool fermi_dns = fermi_bool_dns(sym, dns);
    return {idx_dns_rep, fermi_dns};
  }

  std::pair<idx_t, bool> index_dns_fermi(bit_t dns, int sym, bit_t not_ups,
                                         bit_t fermimask) const {
    bit_t dns_rep = group_action_.apply(sym, dns);
    bit_t dns_rep_c = bitops::extract(dns_rep, not_ups);
    idx_t idx_dns_rep = lintable_dnsc_.index(dns_rep_c);
    bool fermi_dns = (bitops::popcnt(dns & fermimask) & 1);
    fermi_dns ^= fermi_bool_dns(sym, dns);
    return {idx_dns_rep, fermi_dns};
  }

  // index and fermi sign for dns with non-trivial stabilizer
  std::tuple<idx_t, bool, int>
  index_dns_fermi_sym(bit_t dns, gsl::span<int const> syms,
                      gsl::span<bit_t const> dnss_out) const {
    auto [rep_dns, rep_sym] =
        symmetries::representative_sym_subset(dns, group_action_, syms);
    auto it = std::lower_bound(dnss_out.begin(), dnss_out.end(), rep_dns);
    if ((it != dnss_out.end()) && (*it == rep_dns)) {
      bool fermi_dns = fermi_bool_dns(rep_sym, dns);
      return {std::distance(dnss_out.begin(), it), fermi_dns, rep_sym};
    } else {
      return {invalid_index, false, rep_sym};
    }
  }

  std::tuple<idx_t, bool, int>
  index_dns_fermi_sym(bit_t dns, gsl::span<int const> syms,
                      gsl::span<bit_t const> dnss_out, bit_t fermimask) const {
    auto [rep_dns, rep_sym] =
        symmetries::representative_sym_subset(dns, group_action_, syms);
    auto it = std::lower_bound(dnss_out.begin(), dnss_out.end(), rep_dns);
    if ((it != dnss_out.end()) && (*it == rep_dns)) {
      bool fermi_dns = (bitops::popcnt(dns & fermimask) & 1);
      fermi_dns ^= fermi_bool_dns(rep_sym, dns);
      return {std::distance(dnss_out.begin(), it), fermi_dns, rep_sym};
    } else {
      return {invalid_index, false, rep_sym};
    }
  }

  // Retrieving index of representative and symmetries for ups
  idx_t index_ups(bit_t ups) const {
    return idces_up_[lintable_ups_.index(ups)];
  }

  gsl::span<int const> syms_ups(bit_t ups) const {
    idx_t idx_ups = lintable_ups_.index(ups);
    auto [start, length] = sym_limits_up_[idx_ups];
    return {syms_up_.data() + start, length};
  }

  inline std::pair<idx_t, gsl::span<int const>> index_syms_up(bit_t ups) const {
    idx_t idx_ups = lintable_ups_.index(ups);
    idx_t index = idces_up_[idx_ups];
    auto [start, length] = sym_limits_up_[idx_ups];
    return {index, {syms_up_.data() + start, length}};
  }

  // Retrieving dns states and norms for given up configuration
  gsl::span<bit_t const> dns_for_ups_rep(bit_t ups) const {
    idx_t idx_ups = index_ups(ups);
    auto [start, length] = dns_limits_[idx_ups];
    return {dns_storage_.data() + start, length};
  }

  gsl::span<double const> norms_for_ups_rep(bit_t ups) const {
    idx_t idx_ups = index_ups(ups);
    auto [start, length] = dns_limits_[idx_ups];
    return {norms_storage_.data() + start, length};
  }

  std::pair<gsl::span<bit_t const>, gsl::span<double const>>
  dns_norms_for_up_rep(bit_t ups) const {
    idx_t idx_ups = index_ups(ups);
    auto [start, length] = dns_limits_[idx_ups];
    auto dnss = gsl::span<bit_t const>{dns_storage_.data() + start, length};
    auto norms = gsl::span<double const>{norms_storage_.data() + start, length};
    return {dnss, norms};
  }

  // Fermi sign when applying sym on states
  bool fermi_bool_ups(int sym, bit_t ups) const {
    return fermi_bool_ups_table_[sym * raw_ups_size_ +
                                 lintable_ups_.index(ups)];
  }
  bool fermi_bool_dns(int sym, bit_t dns) const {
    return fermi_bool_dns_table_[sym * raw_dns_size_ +
                                 lintable_dns_.index(dns)];
  }

  idx_t dnsc_index(bit_t dns) const { return lintable_dnsc_.index(dns); }

private:
  int n_sites_;
  int n_up_;
  int n_dn_;
  PermutationGroupLookup<bit_t> group_action_;
  Representation irrep_;

  idx_t raw_ups_size_;
  idx_t raw_dns_size_;
  idx_t raw_dnsc_size_;

  indexing::LinTable<bit_t> lintable_ups_;
  indexing::LinTable<bit_t> lintable_dns_;
  indexing::LinTable<bit_t> lintable_dnsc_;

  std::vector<bool> fermi_bool_ups_table_;
  std::vector<bool> fermi_bool_dns_table_;

  std::vector<bit_t> reps_up_;
  std::vector<idx_t> idces_up_;
  std::vector<int> syms_up_;
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits_up_;

  std::vector<idx_t> ups_offset_;
  std::vector<std::pair<span_size_t, span_size_t>> dns_limits_;
  std::vector<bit_t> dns_storage_;
  std::vector<double> norms_storage_;

  idx_t size_;
};

} // namespace hydra::indexing
