// Copyright 2019 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef HYDRA_BASES_BASIS_ELECTRON_H_
#define HYDRA_BASES_BASIS_ELECTRON_H_

#include <hydra/bases/basis_spinhalf.h>
#include <hydra/qns/qn_electron.h>
#include <hydra/states/state_electron.h>

namespace hydra {

template <class bit_t = std_bit_t> class BasisElectronIterator {
public:
  BasisElectronIterator() = default;
  BasisElectronIterator(BasisSpinHalf<bit_t> const &up,
                        BasisSpinHalf<bit_t> const &down,
                        state_electron<bit_t> const &state);

  inline bool operator==(BasisElectronIterator<bit_t> const &rhs) const {
    return ((up_iter_ == rhs.up_iter_) && (down_iter_ == rhs.down_iter_));
  }

  inline bool operator!=(BasisElectronIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
  }

  inline BasisElectronIterator &operator++() {
    ++down_iter_;
    if (down_iter_ == down_end_) {
      down_iter_ = down_begin_;
      ++up_iter_;
    }
    return *this;
  }

  inline state_electron<bit_t> operator*() const {
    return {(*up_iter_).spins, (*down_iter_).spins};
  }

private:
  int n_sites_;
  BasisSpinHalfIterator<bit_t> down_begin_, down_end_;
  BasisSpinHalfIterator<bit_t> down_iter_, up_iter_;
};

/*!
  Class to generate all configurations of a BasisElectron-type Hilbert space
  given the number of upspins and downspins on a given number of sites.

  Usage:
  @code
  #include <hydra/hilberspaces/hubbard.h>

  using BasisElectron = hydra::hilbertspaces::BasisElectron<>;
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::hilbertspaces::Print;

  int n_upspins = 4;
  int n_downspins = 4;
  int n_sites = 8;

  hubbard_qn qn = {n_upspins, n_downspins};
  BasisElectron hs(n_sites, qn);
  int ctr=0;
  for (auto state : hs)
    std::cout << Print(n_sites, state) << std::endl;

  @endcode
  @tparam bit_type elementary datatype to store up and down configurations
                     one of uint16, uint32, uint64
*/
template <class bit_t = std_bit_t> class BasisElectron {
public:
  using qn_t = qn_electron;
  using state_t = state_electron<bit_t>;
  using iterator_t = BasisElectronIterator<bit_t>;

  BasisElectron(int const &n_sites, qn_t const &qn);

  /// returns number of sites of the Hilbert space
  int n_sites() const;

  /// returns the quantum number of the Hilbert space
  qn_t qn() const;

  /// returns iterator to first configuration of Hilbertspace
  iterator_t begin() const;

  /// returns iterator to next-to-last configuration of Hilbertspace
  iterator_t end() const;

  /// returns dimension of the Hilbert space with given quantum numbers
  size_t size() const;

  /// returns raw dimension with all quantum numbers, i.e. pow(4, n_sites)
  size_t rawsize() const;

private:
  number_t n_sites_;
  qn_t qn_;
  iterator_t begin_, end_;
};

} // namespace hydra

#endif