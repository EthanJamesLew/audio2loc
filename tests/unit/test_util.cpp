#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <string>

#include "io.h"
#include "util.h"

using Catch::Approx;
using namespace pl;

TEST_CASE("parse_time reads ISO dates and epoch seconds", "[util]") {
  CHECK(parse_time("1970-01-01T00:00:00Z") == Approx(0.0));
  CHECK(parse_time("2025-09-08") == Approx(1757289600.0));
  CHECK(parse_time("2025-09-08T14:42:00Z") ==
        Approx(1757289600.0 + 14 * 3600 + 42 * 60));
  CHECK(parse_time("2025-09-08T14:42:00.500Z") ==
        Approx(1757289600.0 + 14 * 3600 + 42 * 60 + 0.5));
  CHECK(parse_time("1757289600") == Approx(1757289600.0));
}

TEST_CASE("parse_double rejects trailing garbage", "[util]") {
  double v = 0;
  CHECK(parse_double("3.25", v));
  CHECK(v == 3.25);
  CHECK(parse_double("-1e3", v));
  CHECK(v == -1000);
  CHECK_FALSE(parse_double("12abc", v));
  CHECK_FALSE(parse_double("", v));
}

TEST_CASE("string helpers", "[util]") {
  CHECK(trim("  a b \t\n") == "a b");
  CHECK(lower("PC12") == "pc12");
}

TEST_CASE("the JSON reader handles nesting, escapes and errors", "[io]") {
  JVal j;
  std::string err;
  REQUIRE(parse_json(
      R"({"a":[1,2.5,{"b":"x\"yé"}],"t":true,"n":null,"e":-1e2})", j, &err));
  const JVal* a = j.get("a");
  REQUIRE(a);
  REQUIRE(a->t == JVal::Arr);
  REQUIRE(a->a.size() == 3);
  CHECK(a->a[1].n == 2.5);
  CHECK(a->a[2].get("b")->s == "x\"y\xc3\xa9");
  CHECK(j.get("t")->b);
  CHECK(j.get("n")->t == JVal::Null);
  CHECK(j.get("e")->n == -100);
  CHECK(j.get("missing") == nullptr);
  CHECK_FALSE(parse_json("{\"a\":", j, &err));
  CHECK_FALSE(err.empty());
}

TEST_CASE("json_escape quotes control characters", "[io]") {
  CHECK(json_escape("a\"b\\c") == "a\\\"b\\\\c");
  // Control characters become \u escapes, which any JSON reader accepts.
  JVal j;
  REQUIRE(parse_json("\"" + json_escape("x\ny\t") + "\"", j));
  CHECK(j.s == "x\ny\t");
}
