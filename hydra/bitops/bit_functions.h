#pragma once

#include <limits>

namespace hydra::bitops {
  
template <int bitsize, typename uint_t>
constexpr inline int mask1(((uint_t)1 << bitsize) - 1);

template <typename uint_t>
constexpr inline int nbits(std::numeric_limits<uint_t>::digits);
constexpr inline int floorlog2(int x) noexcept {
  return x == 1 ? 0 : 1 + floorlog2(x >> 1);
}

template <typename uint_t>
constexpr uint_t nbits_mask((floorlog2(nbits<uint_t>) << nbits<uint_t>)-1);
template <typename uint_t>
constexpr uint_t nbits_shift(floorlog2(nbits<uint_t>));

} // namespace hydra::bitops
