#pragma once

#ifdef HYDRA_ENABLE_MPI
#include <mpi.h>
#endif

#include <hydra/operators/bondlist.h>


namespace hydra::testcases::spinhalf {

BondList HBchain(int64_t n_sites, double J1, double J2 = 0);
std::tuple<BondList, arma::vec> HBchain_fullspectrum_nup(int64_t n_sites, int64_t nup);
BondList HB_alltoall(int64_t n_sites);
std::tuple<BondList, double> triangular_12_complex(int64_t nup, double eta);

} // namespace hydra::testcases::spinhalf
