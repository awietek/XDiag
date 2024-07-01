#include <filesystem>
#include <xdiag/all.hpp>

int main(int argc, char **argv) {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;

  say_hello();

  // Parse input arguments
  assert(argc == 7);
  int n_sites = atoi(argv[1]);              // number of sites
  int n_up = atoi(argv[2]);                 // number of upspins
  std::string kname = std::string(argv[3]); // momentum k
  double J = atof(argv[4]);
  double Jd = atof(argv[5]);
  int seed = atoi(argv[6]);

  Log("Diagonalizing H in block nup: {}, k: {}", n_up, kname);

  auto lfile = FileToml(format("checkerboard.{}.JJD.toml", n_sites), 'r');
  std::string odir = format("outfiles/seed.{}", seed);
  std::string ofilename = format(
      "outfile.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.nup.{}.k.{}.seed.{}.h5",
      n_sites, J, Jd, n_up, kname, seed);
  std::filesystem::create_directories(odir);
  auto ofile = FileH5(format("{}/{}", odir, ofilename), "w!");

  auto bonds = BondList(lfile["Interactions"]);
  bonds["J"] = J;
  bonds["JD"] = Jd;
  auto group = PermutationGroup(lfile["Symmetries"]);
  auto irrep = Representation(lfile[kname]);

  Log("Creating block ...");
  tic();
  auto block = Spinhalf(n_sites, n_up, group, irrep);
  toc();
  Log("Dimension: {}", block.size());

  Log("Running Lanczos ...");
  tic();
  int n_eig_to_converge = 3;
  int max_iterations = 300;
  auto tmat = eigvals_lanczos(bonds, block, n_eig_to_converge, 1e-12,
                                  max_iterations, seed, 1e-7);
  toc();

  ofile["Alphas"] = tmat.alphas;
  ofile["Betas"] = tmat.betas;
  ofile["Eigenvalues"] = tmat.eigenvalues;
  ofile["Dimension"] = block.size();

  return EXIT_SUCCESS;
}
