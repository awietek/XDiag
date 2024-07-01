#include <xdiag/all.hpp>

int main(int argc, char **argv) {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;

  say_hello();
  
  // Parse input arguments
  assert(argc == 6);
  int n_sites = atoi(argv[1]);              // number of sites
  int n_up = atoi(argv[2]);                 // number of upspins
  std::string kname = std::string(argv[3]); // momentum k
  double J = atof(argv[4]);
  double Jd = atof(argv[5]);

  Log("Diagonalizing H in block nup: {}, k: {}", n_up, kname);

  auto lfile = FileToml(format("checkerboard.{}.JJD.toml", n_sites), 'r');
  auto ofile = FileH5(
      format(
          "outfiles/outfile.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.nup.{}.k.{}.h5",
          n_sites, J, Jd, n_up, kname),
      "w!");

  auto bonds = BondList(lfile["Interactions"]);
  bonds["J"] = J;
  bonds["JD"] = Jd;
  auto group = PermutationGroup(lfile["Symmetries"]);
  auto irrep_k = Representation(lfile[kname]);

  Log("Creating block ...");
  tic();
  auto block_k = Spinhalf(n_sites, n_up, group, irrep_k);
  toc();
  Log("Dimension: {}", block_k.size());

  Log("Creating matrix ...");
  tic();
  auto H = matrix(bonds, block_k);
  toc();

  Log("Diagonalizing ...");
  tic();
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, H);
  ofile["Eigenvalues"] = eigval;
  toc();

  return EXIT_SUCCESS;
}
