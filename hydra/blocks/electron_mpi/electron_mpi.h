#pragma once

#include <unordered_map>
#include <vector>

#include <hydra/blocks/blocks.h>
#include <hydra/combinatorics/hashes.h>
#include <hydra/common.h>
#include <hydra/indexing/lintable.h>

namespace hydra {

template <class bit_t> class ElectronMPI {

  ElectronMPI() = default;
  ElectronMPI(int n_sites, int nup, int ndn);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }

  inline bool charge_conserved() const { return charge_conserved_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline idx_t size() const { return size_; }
  inline idx_t dim() const { return dim_; }

  inline int process(bit_t spins) const {
    return combinatorics::hash_fnv1(spins) % mpi_size_;
  }

  bool operator==(ElectronMPI const &rhs) const;
  bool operator!=(ElectronMPI const &rhs) const;

private:
  int n_sites_;

  bool charge_conserved_;
  int charge_;
  bool sz_conserved_;
  int sz_;
  int n_up_;
  int n_dn_;

  indexing::LinTable<bit_t> lintable_up_;
  indexing::LinTable<bit_t> lintable_dn_;

  idx_t size_up_;
  idx_t size_dn_;
  idx_t size_;
  idx_t dim_;

  int mpi_rank_;
  int mpi_size_;

  // distribution of spin configurations
  std::vector<bit_t> my_ups_;
  std::unordered_map<bit_t, idx_t> my_ups_offset_;
  std::vector<bit_t> my_dns_;
  std::unordered_map<bit_t, idx_t> my_dns_offset_;
};

} // namespace hydra
