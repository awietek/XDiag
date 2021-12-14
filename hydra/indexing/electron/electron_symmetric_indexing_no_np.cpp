#include "electron_symmetric_indexing_no_np.h"

#include <hydra/blocks/utils/block_utils.h>

#include <hydra/combinatorics/subsets.h>

#include <hydra/symmetries/fermi_bool_table.h>
#include <hydra/symmetries/permutation_group_lookup.h>
#include <hydra/symmetries/representative_list.h>
#include <hydra/symmetries/symmetry_operations.h>

namespace hydra::indexing {

template <class bit_t>
ElectronSymmetricIndexingNoNp<bit_t>::ElectronSymmetricIndexingNoNp(
    int n_sites, PermutationGroup permutation_group, Representation irrep)
    : n_sites_(n_sites),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep),
      raw_ups_size_((idx_t)1 << n_sites), raw_dns_size_((idx_t)1 << n_sites),
      lintable_ups_(n_sites), lintable_dns_(n_sites) {

  using combinatorics::Subsets;

  utils::check_n_sites(n_sites, permutation_group);

  fermi_bool_ups_table_ =
      symmetries::fermi_bool_table(Subsets<bit_t>(n_sites), group_action_);
  fermi_bool_dns_table_ =
      symmetries::fermi_bool_table(Subsets<bit_t>(n_sites), group_action_);

  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          SubsetsIndexing<bit_t>(n_sites), group_action_);

  std::tie(dns_storage_, norms_storage_, dns_limits_, ups_offset_, size_) =
      symmetries::electron_dns_norms_limits_offset_size(
          reps_up_, Subsets<bit_t>(n_sites), group_action_, irrep_);
}

template class ElectronSymmetricIndexingNoNp<uint16_t>;
template class ElectronSymmetricIndexingNoNp<uint32_t>;
template class ElectronSymmetricIndexingNoNp<uint64_t>;

} // namespace hydra::indexing
