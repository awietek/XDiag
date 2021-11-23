#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.h"
#include "../tj/testcases_tj.h"

#include <hydra/all.h>

using namespace hydra;

template <class bit_t>
void test_symmetric_apply(BondList bondlist, Couplings couplings,
                          PermutationGroup space_group,
                          std::vector<Representation> irreps) {
  int n_sites = space_group.n_sites();

  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      if (nup + ndn > n_sites)
        continue;

      for (int k = 0; k < (int)irreps.size(); ++k) {
        auto irrep = irreps[k];
        auto block = tJSymmetric<bit_t>(n_sites, nup, ndn, space_group, irrep);

        if (block.size() > 0) {
          auto H_sym = MatrixCplx(bondlist, couplings, block, block);
          REQUIRE(lila::close(H_sym, lila::Herm(H_sym)));

          // Check whether apply gives the same as matrix multiplication
          auto v = lila::Random<complex>(block.size());
          auto w1 = lila::Mult(H_sym, v);
          auto w2 = lila::ZerosLike(v);
          Apply(bondlist, couplings, block, v, block, w2);
          REQUIRE(lila::close(w1, w2));

          // Compute eigenvalues and compare
          auto evals_mat = lila::EigenvaluesSym(H_sym);
          double e0_mat = evals_mat(0);
          double e0_app = E0Cplx(bondlist, couplings, block);
          // lila::Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
          REQUIRE(std::abs(e0_mat - e0_app) < 1e-7);

          // Compute eigenvalues with real arithmitic
          if (is_real(block.irrep()) && is_real(couplings)) {
            auto H_real = MatrixReal(bondlist, couplings, block, block);
            auto evals_mat_real = lila::EigenvaluesSym(H_real);
            REQUIRE(lila::close(evals_mat_real, evals_mat));

            double e0_mat_real = evals_mat_real(0);
            double e0_app_real = E0Real(bondlist, couplings, block);
            REQUIRE(std::abs(e0_mat_real - e0_app_real) < 1e-7);
          }
        }
      }
    }
  }
}

template <class bit_t> void test_tj_symmetric_apply_chains(int n_sites) {
  using namespace hydra::testcases::tj;
  using namespace hydra::testcases::electron;

  lila::Log.out("Tj chain, symmetric apply test, n_sites: {}", n_sites);
  auto [bondlist, couplings] = tJchain(n_sites, 1.0, 5.0);
  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult(n_sites);
  test_symmetric_apply<uint16_t>(bondlist, couplings, space_group, irreps);
}

TEST_CASE("tj_symmetric_apply", "[blocks]") {
  using namespace hydra::testcases::tj;
  using namespace hydra::testcases::electron;

  // Test linear chains
  for (int n_sites = 2; n_sites < 8; ++n_sites) {
    test_tj_symmetric_apply_chains<uint16_t>(n_sites);
    test_tj_symmetric_apply_chains<uint32_t>(n_sites);
    test_tj_symmetric_apply_chains<uint64_t>(n_sites);
  }
  {
    // test a 3x3 triangular lattice
    lila::Log.out("Tj 3x3 triangular, symmetric apply test");
    std::string lfile = "data/triangular.9.hop.sublattices.tsl.lat";

    auto bondlist = read_bondlist(lfile);
    Couplings couplings;
    couplings["T"] = 1.0;
    couplings["U"] = 5.0;
    auto permutations = hydra::read_permutations(lfile);
    auto space_group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_represenation(lfile, name));
    }
    test_symmetric_apply<uint16_t>(bondlist, couplings, space_group, irreps);
    test_symmetric_apply<uint32_t>(bondlist, couplings, space_group, irreps);
    test_symmetric_apply<uint64_t>(bondlist, couplings, space_group, irreps);
  }

  {
    // test a 3x3 triangular lattice with complex flux
    lila::Log.out("tJ 3x3 triangular staggered, symmetric apply test, complex");
    std::string lfile =
        "data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.lat";

    auto bondlist = read_bondlist(lfile);
    std::vector<double> etas{0.0, 0.1, 0.2, 0.3};
    auto permutations = hydra::read_permutations(lfile);
    auto space_group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    std::vector<int> multiplicities;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_represenation(lfile, name));
      multiplicities.push_back(mult);
    }

    for (auto eta : etas) {
      lila::Log("eta: {:.2f}", eta);
      Couplings couplings;
      couplings["TPHI"] = complex(cos(eta * M_PI), sin(eta * M_PI));
      couplings["JPHI"] = complex(cos(2 * eta * M_PI), sin(2 * eta * M_PI));

      test_symmetric_apply<uint16_t>(bondlist, couplings, space_group, irreps);
      test_symmetric_apply<uint32_t>(bondlist, couplings, space_group, irreps);
      test_symmetric_apply<uint64_t>(bondlist, couplings, space_group, irreps);
    }
  }
}
