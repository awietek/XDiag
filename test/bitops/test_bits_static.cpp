#include <hydra/bitops/bits_static.h>


TEST_CASE("bits_static", "[bitops]") {
  Log.out("Testing BitsStatic");
  Log.set_verbosity(1);
  BitsStatic<2, 1, unsigned char> bits = 12;
  std::cout << bits <<"\n";
}
