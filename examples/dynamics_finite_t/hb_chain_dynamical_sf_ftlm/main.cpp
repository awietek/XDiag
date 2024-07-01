#include <filesystem>
#include <xdiag/all.hpp>

int main(int argc, char **argv) {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;
  using hdf5_opts::append;
  using hdf5_opts::trans;

  // Parse input arguments
  assert(argc == 6);
  int n_sites = atoi(argv[1]); // number of sites
  int n_up = atoi(argv[2]);    // number of upspins
  int k = atoi(argv[3]);       // momentum k
  int seed = atoi(argv[4]);    // random seed
  int n_iters = atoi(argv[5]); // number of iterations

  Log.set_verbosity(1);

  // Define directories for output / scratch data
  std::string outdir = format("outfiles/N.{}/seed.{}", n_sites, seed);
  std::string outfile =
      format("{}/outfile.N.{}.nup.{}.k.{}.seed.{}.iters.{}.h5", outdir, n_sites,
             n_up, k, seed, n_iters);
  std::string scratchdir =
      format("/scratch/martin/tmp/ftlm.N.{}.nup.{}.k.{}.seed.{}.iters.{}/",
             n_sites, n_up, k, seed, n_iters);
  std::filesystem::create_directories(outdir);
  std::filesystem::create_directories(scratchdir);

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
  for (int q = 0; q < n_sites; ++q) {
    complex phase = exp(2i * pi * q / n_sites);
    auto irrep = generated_irrep(perm, phase);
    irreps.push_back(irrep);
  }

  ////////////////////////////////////////////////////////
  Log("Creating Lanczos vector matrix V from |r> ...");
  auto block = Spinhalf(n_sites, n_up, group, irreps[k]);
  ivec dim_k(1);
  dim_k(0) = block.size();
  dim_k.save(hdf5_name(outfile, "dim", append));

  auto rstate = random_state(block, seed);
  auto dump_V = [scratchdir](int iteration, cx_vec const &vec) {
    std::string filename = format("{}/V{}.arm", scratchdir, iteration);
    Log("writing V vector {} -> {}", iteration, filename);
    vec.save(filename);
  };
  auto T = lanczos_vector_apply_inplace(bonds, rstate, dump_V, n_iters);
  T.alphas().save(hdf5_name(outfile, "T_alphas", append));
  T.betas().save(hdf5_name(outfile, "T_betas", append));

  ////////////////////////////////////////////////////////
  Log("Computing ground state ...");
  XDIAG_SHOW(block);
  auto gs = groundstate(bonds, block);

  for (int q = 0; q < n_sites; ++q) {

    Log("Creating S(q)|g.s.> (q={}) ...", q);
    auto S_of_q = symmetrized_operator(Bond("SZ", 0), group, irreps[q]);
    auto block_q = Spinhalf(n_sites, n_up, group, irreps[k] * irreps[q]);
    ivec dim_q(1);
    dim_q(0) = block_q.size();
    dim_q.save(hdf5_name(outfile, format("q_{}/dim", q), append));

    auto w0 = State(block_q);
    apply(S_of_q, gs, w0);

    // Compute norm and store
    vec nrm(1);
    nrm(0) = norm(w0);
    nrm.save(hdf5_name(outfile, format("q_{}/norm", q), append));

    if (nrm(0) < 1e-12) {
      Log("Zero norm of S(q)|g.s.>");
      vec alphas = zeros(0);
      vec betas = zeros(0);
      alphas.save(hdf5_name(outfile, format("q_{}/S_alphas", q), append));
      betas.save(hdf5_name(outfile, format("q_{}/S_betas", q), append));
      continue;
    }

    Log("Creating Lanczos vector matrix W from S(q)|g.s.>...");
    auto dump_W = [scratchdir, q](int iteration, cx_vec const &vec) {
      std::string filename = format("{}/W{}.arm", scratchdir, iteration);
      Log("writing W vector {} (q={}) -> {}", iteration, q, filename);
      vec.save(filename);
    };
    auto S = lanczos_vector_apply_inplace(bonds, w0, dump_W, n_iters);
    S.alphas().save(hdf5_name(outfile, format("q_{}/S_alphas", q), append));
    S.betas().save(hdf5_name(outfile, format("q_{}/S_betas", q), append));

    Log("Computing A = W S(q) V (q={}) ...", q);
    auto A = cx_mat(S.size(), T.size());

    for (int n = 0; n < T.size(); ++n) {
      Log("Computing A (q={}) -> n={}", q, n);
      tic();
      auto v = State(block);
      v.vector().load(format("{}/V{}.arm", scratchdir, n));
      auto S_of_q_v = State(block_q);
      apply(S_of_q, v, S_of_q_v);

      for (int m = 0; m < S.size(); ++m) {
        auto w = State(block_q);
        w.vector().load(format("{}/W{}.arm", scratchdir, m));
        A(m, n) = dot(w, S_of_q_v);
      }
      toc();
    }
    A.save(hdf5_name(outfile, format("q_{}/A", q), append + trans));

    // Remove W vectors
    for (int m = 0; m < S.size(); ++m) {
      std::remove(format("{}/W{}.arm", scratchdir, m).c_str());
    }
  }

  // Remove V vectors
  for (int n = 0; n < T.size(); ++n) {
    std::remove(format("{}/V{}.arm", scratchdir, n).c_str());
  }
  std::filesystem::remove_all(scratchdir.c_str());

  return EXIT_SUCCESS;
}
