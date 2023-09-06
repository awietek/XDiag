#include "bitops.h"

#include <hydra/bitops/bits_static.h>

#include <sstream>

namespace hydra::bitops {

// bits_to_string
template <typename bit_t>
std::string bits_to_string(bit_t bits, int n, bool reverse) {
  std::stringstream s;
  for (int i = 0; i < n; ++i)
    s << gbit(bits, i);
  std::string st = s.str();
  return reverse ? std::string(st.rbegin(), st.rend()) : st;
}

template std::string bits_to_string(int16_t, int, bool);
template std::string bits_to_string(int32_t, int, bool);
template std::string bits_to_string(int64_t, int, bool);

template std::string bits_to_string(uint16_t, int, bool);
template std::string bits_to_string(uint32_t, int, bool);
template std::string bits_to_string(uint64_t, int, bool);

template std::string bits_to_string(BitsStatic<1, 1, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<1, 1, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<1, 1, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<1, 2, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<1, 2, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<1, 2, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<1, 3, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<1, 3, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<1, 3, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<1, 4, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<1, 4, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<1, 4, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<2, 1, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<2, 1, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<2, 1, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<2, 2, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<2, 2, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<2, 2, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<2, 3, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<2, 3, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<2, 3, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<2, 4, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<2, 4, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<2, 4, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<4, 1, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<4, 1, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<4, 1, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<4, 2, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<4, 2, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<4, 2, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<4, 3, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<4, 3, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<4, 3, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<4, 4, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<4, 4, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<4, 4, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<8, 1, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<8, 1, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<8, 1, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<8, 2, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<8, 2, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<8, 2, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<8, 3, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<8, 3, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<8, 3, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<8, 4, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<8, 4, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<8, 4, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<16, 1, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<16, 1, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<16, 1, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<16, 2, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<16, 2, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<16, 2, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<16, 3, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<16, 3, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<16, 3, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<16, 4, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<16, 4, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<16, 4, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<32, 1, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<32, 1, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<32, 1, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<32, 2, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<32, 2, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<32, 2, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<32, 3, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<32, 3, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<32, 3, uint64_t>, int, bool);

template std::string bits_to_string(BitsStatic<32, 4, uint16_t>, int, bool);
template std::string bits_to_string(BitsStatic<32, 4, uint32_t>, int, bool);
template std::string bits_to_string(BitsStatic<32, 4, uint64_t>, int, bool);

} // namespace hydra::bitops
