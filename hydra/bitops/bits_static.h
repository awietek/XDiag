#pragma once
#include <limits>

#include <hydra/common.h>

namespace hydra::bitops {


// three template parameters:
// nwords: number of uints used to code the state
// bitsize: numbers of physical bits used for a logical bit
// uint_t: type of the underlying unsigned integer
template <int nwords, int bitsize, typename uint_t> class BitsStatic {
public:
  BitsStatic() = default;
  inline BitsStatic(int n) { words_[0] = (uint_t)n; }
  inline BitsStatic(uint_t n) { words_[0] = n; }
  inline uint_t word(int n) const { return words_[n]; }

  inline BitsStatic operator<<(int shift) {
    int nperm = (shift * bitsize) / nbits<uint_t>;
    int nshift = (shift * bitsize) % nbits<uint_t>;
    int nshift_anti = nbits<uint_t> - nshift;

    BitsStatic res;
    for (int i = 0; i < nwords - 1 - nperm; ++i) {
      res.words_[i + nperm] |= words_[i] << nshift;
      res.words_[i + 1 + nperm] |= words_[i] >> nshift_anti;
    }
  }

private:
  std::array<uint_t, nwords> words_;
};


} // namespace hydra::bitops
