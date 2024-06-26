#include "../catch.hpp"

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/io/args.hpp>

TEST_CASE("args", "[io]") {
  using namespace xdiag;
  using namespace arma;
  Log("testing Args");

  Args args({{"bool", true},
             {"string", "hello"},
             {"int", 42},
             {"double", 1.23},
             {"complex", (complex)(1.23 + 3.21i)}});

  REQUIRE(args.defined("bool"));
  REQUIRE(args.defined("string"));
  REQUIRE(args.defined("int"));
  REQUIRE(args.defined("double"));
  REQUIRE(args.defined("complex"));
  REQUIRE(!args.defined("vec"));

  REQUIRE(args["bool"].as<bool>() == true);
  REQUIRE(args["string"].as<std::string>() == "hello");
  REQUIRE(args["int"].as<int>() == 42);
  REQUIRE(args["double"].as<double>() == 1.23);
  REQUIRE(args["complex"].as<complex>() == (complex)(1.23 + 3.21i));

  REQUIRE(args["asdf0"].as<bool>(false) == false);
  REQUIRE(args["asdf1"].as<std::string>("jjj") == "jjj");
  REQUIRE(args["asdf2"].as<int>(43) == 43);
  REQUIRE(args["asdf3"].as<double>(4.3) == 4.3);
  REQUIRE(args["asdf4"].as<complex>(4.3 + 1.2i) == (complex)(4.3 + 1.2i));
}
