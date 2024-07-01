#include <filesystem>
#include <xdiag/all.hpp>

int main(int argc, char** argv) {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;
  using hdf5_opts::append;
  using hdf5_opts::trans;

  // Parse input arguments
  assert(argc == 4);
  int n_sites = atoi(argv[1]); // number of sites
  int n_up = atoi(argv[2]);    // number of upspins
  int k = atoi(argv[3]);       // momentum k

  say_hello();
  
  // Define directory / file to store output data
  std::string outdir = format("outfiles/N.{}", n_sites);
  std::string outfile = format("{}/outfile.N.{}.nup.{}.k.{}.h5",
			       outdir, n_sites, n_up, k);
  std::filesystem::create_directories(outdir);

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

  // Create the irrep at momentum k
  complex phase = exp(2i * pi * k / (double)n_sites);
  auto irrep = generated_irrep(perm, phase);
  
  auto block = Spinhalf(n_sites, n_up, group, irrep);
  XDIAG_SHOW(block);
  
  // Compute eigendecomposition of Hamiltonian
  Log("Creating H");
  cx_mat H = matrixC(bonds, block);

  Log("Diagonalizing H");
  vec eigval;
  cx_mat eigvec;
  eig_sym(eigval, eigvec, H);
  eigval.save(hdf5_name(outfile, "eigenvalues", append));

  // Loop over different momenta q
  for (int q = 0; q <= n_sites; ++q) {

    Log("Computing <n|S(q)|m> (q={})", q);
    
    // Create S(q) operator
    BondList S_of_q_bonds;
    for (int s = 0; s < n_sites; ++s) {
      complex phase = exp(2i * pi * q * s / (double)n_sites);
      S_of_q_bonds << Bond("SZ", phase / n_sites, s);
    }

    // Create block at momentum k + q
    complex phase_q = exp(2i * pi * (k+q) / (double)n_sites);
    auto irrep_q = generated_irrep(perm, phase_q);
    auto block_q = Spinhalf(n_sites, n_up, group, irrep_q);
    cx_mat H_q = matrixC(bonds, block_q);

    vec eigval_q;
    cx_mat eigvec_q;
    eig_sym(eigval_q, eigvec_q, H_q);
    eigval.save(hdf5_name(outfile, "eigenvalues", append));
    
    // Compute matrix elements of S(q)
    cx_mat S_of_q = matrixC(S_of_q_bonds, block, block_q);
    cx_mat S_of_q_eig = eigvec_q.t() * S_of_q * eigvec;
    S_of_q_eig.save(hdf5_name(outfile, format("S_of_q_{}_eig", q), append + trans));
  }

  return EXIT_SUCCESS;
}
