#include "indexing_symmetric_fermionic.hpp"

#include <xdiag/blocks/utils/block_utils.hpp>
#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/symmetries/fermi_sign.hpp>
#include <xdiag/symmetries/permutation_group_action.hpp>
#include <xdiag/symmetries/permutation_group_lookup.hpp>
#include <xdiag/symmetries/symmetry_operations.hpp>

namespace xdiag::indexing {

template <typename bit_t, class GroupAction>
IndexingSymmetricFermionic<bit_t, GroupAction>::IndexingSymmetricFermionic(
    int n_sites, int n_fermions, PermutationGroup permutation_group,
    Representation irrep)
    : lin_table_(n_sites, n_fermions),
      index_of_raw_state_(combinatorics::binomial(n_sites, n_fermions),
                          invalid_index),
      norm_of_raw_state_(combinatorics::binomial(n_sites, n_fermions), 0.) {
  using combinatorics::Combinations;

  utils::check_nup_spinhalf(n_sites, n_fermions, "IndexingSymmetricFermionic");

  // if not all symmetries are allowed by irrep, choose a subgroup
  if (irrep.allowed_symmetries().size() > 0) {
    permutation_group = permutation_group.subgroup(irrep.allowed_symmetries());
  }
  auto group_action = GroupAction(permutation_group);

  // Go through non symmetrized states and register representatives
  int64_t idx = 0;
  int64_t n_representatives = 0;
  std::vector<int> fermi_work(n_sites, 0);

  for (bit_t state : Combinations(n_sites, n_fermions)) {

    bit_t rep = symmetries::representative(state, group_action);
    // register state if it's a representative
    if (rep == state) {

      double norm = symmetries::norm_fermionic(rep, group_action, irrep);

      if (norm > 1e-6) { // tolerance big as 1e-6 since root is taken
        int64_t idx = lin_table_.index(rep);
        states_.push_back(rep);
        norm_of_raw_state_[idx] = norm;
        index_of_raw_state_[idx] = n_representatives++;
      }
    }
    ++idx;
  }

  // Go through non symmetrized states and fill indices for all states
  idx = 0;
  const int *sym_ptr = permutation_group.permutation_array().data();
  for (bit_t state : Combinations(n_sites, n_fermions)) {

    // Compute the representative of state and the corresponding symmetry
    auto [rep, rep_sym] = symmetries::representative_sym(state, group_action);

    if (rep != state) {
      complex norm = norm_of_raw_state_[lin_table_.index(rep)];
      if (std::abs(norm) > 1e-6) {
        norm_of_raw_state_[idx] =
            symmetries::fermi_sign_of_permutation(
                state, sym_ptr + rep_sym * n_sites, fermi_work.data()) *
            irrep.character(rep_sym) * norm;
        index_of_raw_state_[idx] = index_of_raw_state_[lin_table_.index(rep)];
      }
    }
    ++idx;
  }

  size_ = (int64_t)states_.size();
}

template class IndexingSymmetricFermionic<uint16_t, PermutationGroupAction>;
template class IndexingSymmetricFermionic<uint32_t, PermutationGroupAction>;
template class IndexingSymmetricFermionic<uint64_t, PermutationGroupAction>;

template class IndexingSymmetricFermionic<uint16_t,
                                          PermutationGroupLookup<uint16_t>>;
template class IndexingSymmetricFermionic<uint32_t,
                                          PermutationGroupLookup<uint32_t>>;
template class IndexingSymmetricFermionic<uint64_t,
                                          PermutationGroupLookup<uint64_t>>;

} // namespace xdiag::indexing
