#include "../../catch.hpp"

#include <iostream>

#include <hydra/all.h>

using namespace hydra;
using namespace hydra::combinatorics;

template <class bit_t>
void test_indices(ElectronSymmetric<bit_t> const &electron) {
  // int n_sites = electron.n_sites();

  idx_t idx = 0;
  for (bit_t ups : electron.reps_up_) {
    auto const &dnss = electron.dns_for_up_rep(ups);
    for (bit_t dns : dnss) {
      // lila::Log.out("{} {}", bits_to_string(ups, n_sites),
      // bits_to_string(dns, n_sites));
      idx_t idx2 = electron.index(ups, dns);
      REQUIRE(idx == idx2);
      ++idx;
    }
  }
}

template <class bit_t> void test_electron_chain(int n_sites) {
  using bitops::bits_to_string;

  lila::Log("electron_symmetric: chain block test. N: {}", n_sites);

  // cyclic group as space group
  std::vector<std::vector<int>> permutations;
  for (int sym = 0; sym < n_sites; ++sym) {
    std::vector<int> permutation;
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation.push_back(newsite);
    }
    permutations.push_back(permutation);
  }
  auto space_group = PermutationGroup(permutations);

  // Create irrep with K momentum
  idx_t sum_of_dims = 0;

  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      idx_t sum_of_dims_updn = 0;

      for (int k = 0; k < n_sites; ++k) {

        // lila::Log.out("N: {}, nup: {}, ndn: {}, k:{} ", n_sites, nup, ndn,
        // k);

        // Create characters
        std::vector<complex> chis;
        for (int l = 0; l < n_sites; ++l)
          chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                          std::sin(2 * M_PI * l * k / n_sites)});
        auto irrep = Representation(chis);

        auto electron2 =
            ElectronSymmetric<bit_t>(n_sites, nup, ndn, space_group, irrep);

        sum_of_dims += electron2.size();
        sum_of_dims_updn += electron2.size();
      }
      CHECK(sum_of_dims_updn ==
            binomial(n_sites, nup) * binomial(n_sites, ndn));
    }
  }
  REQUIRE(sum_of_dims == pow(4, n_sites));
}

TEST_CASE("electron_symmetric", "[blocks][electron_symmetric]") {

  // Test the Hubbard chain
  for (int n_sites = 1; n_sites < 7; ++n_sites) {
    test_electron_chain<uint16_t>(n_sites);
    test_electron_chain<uint32_t>(n_sites);
    test_electron_chain<uint64_t>(n_sites);
  }

  // test a 3x3 triangular lattice
  int n_sites = 9;
  lila::Log.out("electron_symmetric: triangular 3x3 block test");

  std::vector<std::pair<std::string, int>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};

  std::string lfile = "data/triangular.9.hop.sublattices.tsl.lat";
  auto permutations = hydra::read_permutations(lfile);
  auto space_group = PermutationGroup(permutations);

  idx_t sum_dim = 0;
  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      idx_t sum_dim_updn = 0;

      for (auto [name, mult] : rep_name_mult) {
        auto irrep = read_represenation(lfile, name);
        auto electron =
            ElectronSymmetric(n_sites, nup, ndn, space_group, irrep);

        idx_t dim = electron.size() * mult;
        // lila::Log.out(
        //     "Hubbard Triangular 3x3: n_sites: {}, nup: {}, ndn: {}, k: "
        //     "{}(x{}), size: {}",
        //     n_sites, nup, ndn, name, mult, electron.size());
        test_indices(electron);

        sum_dim_updn += dim;
        sum_dim += dim;
      }
      auto electron_nosym = Electron(n_sites, nup, ndn);
      REQUIRE(sum_dim_updn == electron_nosym.size());
      REQUIRE(sum_dim_updn == binomial(n_sites, nup) * binomial(n_sites, ndn));

      // lila::Log.out("size: {} {}", sum_dim_updn, electron_nosym.size());
    }
  }
  REQUIRE(sum_dim == (idx_t)pow(4, n_sites));
}
