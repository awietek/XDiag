#include <hydra/all.h>

int main() {
  using namespace hydra;

  // Define the number of lattice sites, number of upspins
  int n_sites = 8;
  int n_up = 4;

  // Create the corresponding block of the Spinhalf Hilbertspace 
  auto block = Spinhalf(n_sites, n_up);

  // Define Heisenberg chain with interaction-type HB and coupling J
  BondList bonds;
  for (int s = 0; s < n_sites; ++s) {
    bonds << Bond("HB", "J", {s, (s + 1) % n_sites});
  }

  // Set the coupling J to 1.0
  Couplings couplings;
  couplings["J"] = 1.0;

  // Compute the full Hamiltonian matrix (real, double precision)
  auto H = MatrixReal(bonds, couplings, block, block);

  // Perform a full diagonalization
  auto eigs = lila::EigenvaluesSym(H);

  // Output the eigenvalues
  LilaPrint(eigs);
  
  return EXIT_SUCCESS;
}
