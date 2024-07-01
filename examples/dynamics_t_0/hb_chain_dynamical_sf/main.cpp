#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;

  say_hello();
  
  int n_sites = 16;
  int n_up = n_sites / 2;

  // Create nearest-neighbor Heisenberg model
  BondList bonds;
  for (int s = 0; s < n_sites; ++s) {
    bonds << Bond("HB", "J", {s, (s + 1) % n_sites});
  }
  bonds["J"] = 1.0;

  // Create the permutation group
  std::vector<int> translation;
  for (int s = 0; s < n_sites; ++s) {
    translation.push_back((s + 1) % n_sites);
  }
  Permutation perm(translation);
  auto group = generated_group(perm);

  // Create the irreps at momenta k
  std::vector<Representation> irreps;
  for (int k = 0; k < n_sites; ++k) {
    complex phase = exp(2i * pi * k / n_sites);
    auto irrep = generated_irrep(perm, phase);
    irreps.push_back(irrep);
  }

  // compute groundstate (known to be at k=0)
  Log("Computing ground state ...");
  auto block = Spinhalf(n_sites, n_up, group, irreps[0]);
  auto [e0, gs] = eig0(bonds, block);
  // auto gs = resgs.eigenvectors.col(0);
  gs.make_complex();
  Log("done.");

  for (int q = 0; q < n_sites; ++q) {

    Log("Dynamical Lanczos iterations for q={}", q);

    // Compute S(q) |g.s.>
    auto S_of_q = symmetrized_operator(Bond("SZ", 0), group, irreps[q]);
    auto block_q = Spinhalf(n_sites, n_up, group, irreps[q]);
    auto v0 = State(block_q);
    v0.make_complex();
    apply(S_of_q, gs, v0);


    double nrm = norm(v0);
    v0 /= nrm;

    // Perform 200 Lanczos iterations for dynamics starting from v0
    auto res = eigvals_lanczos(bonds, block_q, v0, 1, 0., 200, true, 1e-7);
    auto alphas = res.alphas;
    auto betas = res.betas;

    // Write alphas, betas, and norm to file for further processing
    alphas.save(
        fmt::format("outfiles/alphas.N.{}.nup.{}.q.{}.txt", n_sites, n_up, q),
        raw_ascii);
    betas.save(
        fmt::format("outfiles/betas.N.{}.nup.{}.q.{}.txt", n_sites, n_up, q),
        raw_ascii);
    std::ofstream of;
    of.open(
        fmt::format("outfiles/norm.N.{}.nup.{}.q.{}.txt", n_sites, n_up, q));
    of << nrm;
    of.close();
  }

  return EXIT_SUCCESS;
}
