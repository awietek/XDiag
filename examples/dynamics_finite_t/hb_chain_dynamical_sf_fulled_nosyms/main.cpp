#include <filesystem>
#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;
  using hdf5_opts::append;

  say_hello();
  
  // Parse input arguments
  int n_sites = 12;

  // Define directory / file to store output data
  std::string outdir = format("outfiles/N.{}", n_sites);
  std::string outfile = format("{}/outfile.N.{}.h5", outdir, n_sites);
  std::filesystem::create_directories(outdir);

  // Create nearest-neighbor Heisenberg model
  BondList bonds;
  for (int s = 0; s < n_sites; ++s) {
    bonds << Bond("HB", "J", {s, (s + 1) % n_sites});
  }
  bonds["J"] = 1.0;

  auto block = Spinhalf(n_sites);

  // Compute eigendecomposition of Hamiltonian
  Log("Creating H");
  mat H = matrix(bonds, block);

  Log("Diagonalizing H");
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, H);
  eigval.save(hdf5_name(outfile, "eigenvalues", append));

  // Loop over different momenta q
  for (int q = 0; q <= n_sites; ++q) {

    Log("Computing <n|S(q)|m> (q={})", q);
    
    // Create S(q) operator
    BondList S_of_q_bonds;
    for (int s = 0; s < n_sites; ++s) {
      complex phase = exp(2i * pi * q * s / n_sites);
      S_of_q_bonds << Bond("SZ", phase / n_sites, s);
    }

    // Compute matrix elements of S(q)
    cx_mat S_of_q = matrixC(S_of_q_bonds, block);
    cx_mat S_of_q_eig = eigvec.t() * S_of_q * eigvec;
    S_of_q_eig.save(hdf5_name(outfile, format("S_of_q_{}_eig", q), append));
  }

  return EXIT_SUCCESS;
}
