#include "../../catch.hpp"

#include <iostream>

#include "testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;

TEST_CASE("electron_apply", "[blocks][electron]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;
  Couplings couplings;

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {

    lila::Log("electron_apply: Hubbard random all-to-all test (real), N: {}",
              n_sites);
    std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);
    couplings["U"] = 5.0;

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron<uint32_t>(n_sites, nup, ndn);
        auto H = MatrixReal(bondlist, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        // Check whether apply gives the same as matrix multiplication
        auto v = lila::Random<double>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        Apply(bondlist, couplings, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        // Compute eigenvalues and compare
        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = E0Real(bondlist, couplings, block);
        // lila::Log.out("nup: {}, ndn: {}, e0_mat: {}, e0_app: {}", nup, ndn,
        // e0_mat, e0_app);
        CHECK(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
      }
  }

  // /////////////////
  // Test Fermion all to all, free fermions (cplx, up/dn different)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    lila::Log("electron_apply: Hubbard random all-to-all test (cplx), N: {}",
              n_sites);
    std::tie(bondlist, couplings) = freefermion_alltoall_complex_updn(n_sites);
    couplings["U"] = 5.0;

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron<uint32_t>(n_sites, nup, ndn);
        auto H = MatrixCplx(bondlist, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        // Check whether apply gives the same as matrix multiplication
        auto v = lila::Random<complex>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        Apply(bondlist, couplings, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        // Compute eigenvalues and compare
        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = E0Cplx(bondlist, couplings, block);
        // lila::Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
      }
  }

  /////////////////////////////////////////////////////
  // Henry's MATLAB code test (tests Heisenberg terms)
  lila::Log("electron_apply: U-hopping-HB apply of Henry's Matlab code");
  {
    auto [bondlist, couplings, eigs_correct] = randomAlltoAll4NoU();

    int n_sites = 4;
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {
        auto block = Electron(n_sites, nup, ndn);
        auto H = MatrixReal(bondlist, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        // Check whether apply gives the same as matrix multiplication
        auto v = lila::Random<double>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        Apply(bondlist, couplings, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        // Compute eigenvalues and compare
        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = E0Cplx(bondlist, couplings, block);
        // lila::Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
      }

    std::tie(bondlist, couplings, eigs_correct) = randomAlltoAll4();
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {
        auto block = Electron(n_sites, nup, ndn);
        auto H = MatrixReal(bondlist, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        // Check whether apply gives the same as matrix multiplication
        auto v = lila::Random<double>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        Apply(bondlist, couplings, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        // Compute eigenvalues and compare
        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = E0Cplx(bondlist, couplings, block);
        // lila::Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
      }
  }
}
