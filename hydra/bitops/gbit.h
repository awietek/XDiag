#pragma once

#include <cstdint>

// (un)comment to use intrinsic BMI1 bextr instruction
// #define HAS_BEXTR

#if defined(__BMI__) && defined(HAS_BEXTR)
#include <immintrin.h>
#endif

namespace hydra::bitops {

#if defined(__BMI__) && defined(HAS_BEXTR)
constexpr uint32_t bextr(uint16_t src, uint32_t start, uint32_t len) noexcept {
  return _bextr_u32((uint32_t)src, start, len);
}
constexpr uint32_t bextr(uint32_t src, uint32_t start, uint32_t len) noexcept {
  return _bextr_u32(src, start, len);
}
constexpr uint64_t bextr(uint64_t src, uint32_t start, uint32_t len) noexcept {
  return _bextr_u64(src, start, len);
}
#endif


// gbit
template <class bit_t> constexpr bit_t gbit(bit_t x, uint32_t n) noexcept {
#if defined(__BMI__) && defined(HAS_BEXTR)
  return bextr(x, n, 1);
#else
  return (x >> n) & (bit_t)1;
#endif
}

// gbits
template <class bit_t>
constexpr bit_t gbits(bit_t src, uint32_t start, uint32_t len) noexcept {
#if defined(__BMI__) && defined(HAS_BEXTR)
  return bextr(src, start, len);
#else
  return (src >> start) & (((bit_t)1 << len) - 1);
#endif
}

} // namespace hydra::bitops
