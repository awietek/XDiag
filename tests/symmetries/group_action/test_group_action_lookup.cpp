#include "../../catch.hpp"

#include <iostream>

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>
#include <xdiag/symmetries/group_action/group_action_lookup.hpp>

using namespace xdiag;

template <class bit_t> void test_permutation_group_lookup(int64_t n_sites) {
  using combinatorics::Subsets;

  // test cyclic group
  std::vector<Permutation> permutation_array;
  for (int64_t sym = 0; sym < n_sites; ++sym) {

    std::vector<int64_t> pv;
    for (int64_t site = 0; site < n_sites; ++site) {
      int64_t newsite = (site + sym) % n_sites;
      pv.push_back(newsite);
    }
    permutation_array.push_back(Permutation(pv));
  }
  auto perm_group = PermutationGroup(permutation_array);
  auto action = GroupAction(perm_group);
  auto lookup = GroupActionLookup<bit_t>(perm_group);

  // Check whether representative is smallest in orbit
  for (int64_t sym = 0; sym < perm_group.size(); ++sym)
    for (bit_t state : Subsets<bit_t>(n_sites)) {

      auto s1 = action.apply(sym, state);
      auto s2 = lookup.apply(sym, state);
      CHECK(s1 == s2);
    }
}

TEST_CASE("GroupActionLookup", "[symmetries]") {
  xdiag::Log.out("PermutationGroupLookup <-> GroupAction comparison");
  for (int64_t n_sites = 1; n_sites < 6; ++n_sites) {
    test_permutation_group_lookup<uint16_t>(n_sites);
    test_permutation_group_lookup<uint32_t>(n_sites);
    test_permutation_group_lookup<uint64_t>(n_sites);
  }
  xdiag::Log("done");
}
