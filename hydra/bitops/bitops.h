#pragma once

#include <hydra/bitops/bit_functions.h>
#include <hydra/bitops/extract_deposit.h>
#include <hydra/bitops/gbit.h>
#include <hydra/bitops/popcnt.h>

#include <string>

// #if defined(__BMI__)
// #warning has BMI1
// #endif

// #if defined(__BMI2__)
// #warning has BMI2
// #endif

namespace hydra::bitops {
  
// bits_to_string
template <typename bit_t>
std::string bits_to_string(bit_t bits, int n, bool reverse = true);

} // namespace hydra::bitops
