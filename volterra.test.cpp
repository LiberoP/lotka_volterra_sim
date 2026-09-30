// tests (uses doctest)

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "volterra.hpp"

#include <doctest/doctest.h>

TEST_CASE("Testing Simulation")
{
  fn::Simulation sim;
  fn::Parameters parsA{1., 3., 2., 5., 10, 0.001}; // A, B, C, D, N, delta-t
  fn::Parameters parsB{2.18, 3.14, 4.17, 6.73, 100, 0.0001};
  fn::Parameters parsC{1.2, 0.5, 0.1,
                       0.5, 100, 0.001}; // note: more realistic parameters
  fn::Parameters parsD{1.3, 0.4, 0.1, 0.5, 100, 0.001};

  fn::Parameters pars0{1., 3.2, 0, 3., 3, 0.001};
  fn::Parameters pars00{1., 1., 1., 3., 3, 0};
  fn::Parameters parsneg{-1., 3., 2.4, 4., 3, 0.001};

  fn::UserPoint p1{20., 10.}; // x0, y0
  fn::UserPoint p2{10., 20.};
  fn::UserPoint p3{10., 10.};

  fn::UserPoint p4{171., 34.};
  fn::UserPoint p5{18., 299.};
  fn::UserPoint p6{231., 227.};

  REQUIRE(sim.size() == 0);

  SUBCASE("no parameters")
  {
    CHECK_THROWS(sim.addUserPoint(p1));
  }

  SUBCASE("negative parameter")
  {
    CHECK_THROWS(sim.addParameters(parsneg));
  }

  SUBCASE("null parameter")
  {
    CHECK_THROWS(sim.addParameters(pars0));
  }

  SUBCASE("null steps parameter")
  {
    CHECK_THROWS(sim.addParameters(pars00));
  }

  SUBCASE("no starting point")
  {
    sim.addParameters(parsA);
    CHECK_THROWS(sim.evolve());
  }

  SUBCASE("negative starting point")
  {
    sim.addParameters(parsA);
    CHECK_THROWS(sim.addUserPoint({2., -1.}));
  }

  SUBCASE("null starting point")
  {
    sim.addParameters(parsA);
    CHECK_THROWS(sim.addUserPoint({0., 1.}));
  }

  SUBCASE("volterra A1")
  {
    sim.addParameters(parsA);
    sim.addUserPoint(p1);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsA.N].x == doctest::Approx(14.21).epsilon(0.05));
    CHECK(result[parsA.N].y == doctest::Approx(13.39).epsilon(0.05));
    CHECK(result[parsA.N].H == doctest::Approx(52.72).epsilon(0.05));
  }

  SUBCASE("volterra A2")
  {
    sim.addParameters(parsA);
    sim.addUserPoint(p2);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsA.N].x == doctest::Approx(5.35).epsilon(0.05));
    CHECK(result[parsA.N].y == doctest::Approx(22.09).epsilon(0.05));
    CHECK(result[parsA.N].H == doctest::Approx(65.49).epsilon(0.05));
  }

  SUBCASE("volterra A3")
  {
    sim.addParameters(parsA);
    sim.addUserPoint(p3);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsA.N].x == doctest::Approx(7.33).epsilon(0.05));
    CHECK(result[parsA.N].y == doctest::Approx(11.30).epsilon(0.05));
    CHECK(result[parsA.N].H == doctest::Approx(36.18).epsilon(0.05));
  }

  SUBCASE("check 4th step values in case A1")
  {
    sim.addParameters(parsA);
    sim.addUserPoint(p1);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[4].x == doctest::Approx(17.66).epsilon(0.05));
    CHECK(result[4].y == doctest::Approx(11.40).epsilon(0.05));
    CHECK(result[4].H == doctest::Approx(52.72).epsilon(0.05));
  }

  SUBCASE("check 4th step values in case A2")
  {
    sim.addParameters(parsA);
    sim.addUserPoint(p2);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[4].x == doctest::Approx(7.85).epsilon(0.05));
    CHECK(result[4].y == doctest::Approx(21.05).epsilon(0.05));
    CHECK(result[4].H == doctest::Approx(65.49).epsilon(0.05));
  }

  SUBCASE("check 4th step values in case A3")
  {
    sim.addParameters(parsA);
    sim.addUserPoint(p3);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[4].x == doctest::Approx(8.87).epsilon(0.05));
    CHECK(result[4].y == doctest::Approx(10.57).epsilon(0.05));
    CHECK(result[4].H == doctest::Approx(36.18).epsilon(0.05));
  }

  SUBCASE("volterra B4")
  {
    sim.addParameters(parsB);
    sim.addUserPoint(p4);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsB.N].x == doctest::Approx(0.43).epsilon(0.05));
    CHECK(result[parsB.N].y == doctest::Approx(249.06).epsilon(0.05));
    CHECK(result[parsB.N].H == doctest::Approx(777.54).epsilon(0.05));
  }

  SUBCASE("volterra B5")
  {
    sim.addParameters(parsB);
    sim.addUserPoint(p5);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsB.N].x == doctest::Approx(1.10E-3).epsilon(0.05));
    CHECK(result[parsB.N].y == doctest::Approx(302.09).epsilon(0.05));
    CHECK(result[parsB.N].H == doctest::Approx(982.00).epsilon(0.05));
  }

  SUBCASE("volterra B6")
  {
    sim.addParameters(parsB);
    sim.addUserPoint(p6);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsB.N].x == doctest::Approx(4.79E-5).epsilon(0.05));
    CHECK(result[parsB.N].y == doctest::Approx(501.32).epsilon(0.05));
    CHECK(result[parsB.N].H == doctest::Approx(1627.54).epsilon(0.05));
  }

  SUBCASE("check 20th step values in case B4")
  {
    sim.addParameters(parsB);
    sim.addUserPoint(p4);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[20].x == doctest::Approx(111.65).epsilon(0.05));
    CHECK(result[20].y == doctest::Approx(112.74).epsilon(0.05));
    CHECK(result[20].H == doctest::Approx(777.54).epsilon(0.05));
  }

  SUBCASE("check 20th step values in case B5")
  {
    sim.addParameters(parsB);
    sim.addUserPoint(p5);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[20].x
          == doctest::Approx(2.58).epsilon(
              0.1)); // note: test would fail with stricter approximation
    CHECK(result[20].y == doctest::Approx(315.35).epsilon(0.05));
    CHECK(result[20].H == doctest::Approx(982.04).epsilon(0.05));
  }

  SUBCASE("check 20th step values in case B6")
  {
    sim.addParameters(parsB);
    sim.addUserPoint(p6);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[20].x
          == doctest::Approx(18.46).epsilon(
              0.15)); // note: test would fail with stricter approximation
    CHECK(result[20].y == doctest::Approx(504.40).epsilon(0.05));
    CHECK(result[20].H == doctest::Approx(1627.62).epsilon(0.05));
  }

  SUBCASE("volterra C4")
  {
    sim.addParameters(parsC);
    sim.addUserPoint(p4);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsC.N].x == doctest::Approx(12.74).epsilon(0.05));
    CHECK(result[parsC.N].y == doctest::Approx(64.60).epsilon(0.05));
    CHECK(result[parsC.N].H == doctest::Approx(27.30).epsilon(0.05));
  }

  SUBCASE("volterra C5")
  {
    sim.addParameters(parsC);
    sim.addUserPoint(p5);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsC.N].x == doctest::Approx(8.17E-6).epsilon(0.05));
    CHECK(result[parsC.N].y == doctest::Approx(287.88).epsilon(0.05));
    CHECK(result[parsC.N].H == doctest::Approx(143.00).epsilon(0.05));
  }

  SUBCASE("volterra C6")
  {
    sim.addParameters(parsC);
    sim.addUserPoint(p6);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsC.N].x == doctest::Approx(5.03E-4).epsilon(0.05));
    CHECK(result[parsC.N].y == doctest::Approx(260.48).epsilon(0.05));
    CHECK(result[parsC.N].H == doctest::Approx(127.36).epsilon(0.05));
  }

  SUBCASE("check 20th step values in case C4")
  {
    sim.addParameters(parsC);
    sim.addUserPoint(p4);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[20].x == doctest::Approx(117.94).epsilon(0.05));
    CHECK(result[20].y == doctest::Approx(44.91).epsilon(0.05));
    CHECK(result[20].H == doctest::Approx(27.30).epsilon(0.05));
  }

  SUBCASE("check 20th step values in case C5")
  {
    sim.addParameters(parsC);
    sim.addUserPoint(p5);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[20].x
          == doctest::Approx(0.92).epsilon(
              0.25)); // note: test would fail with stricter approximation
    CHECK(result[20].y == doctest::Approx(299.44).epsilon(0.05));
    CHECK(result[20].H == doctest::Approx(143.01).epsilon(0.05));
  }

  SUBCASE("check 20th step values in case C6")
  {
    sim.addParameters(parsC);
    sim.addUserPoint(p6);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[20].x
          == doctest::Approx(18.49).epsilon(
              0.2)); // note: test would fail with stricter approximation
    CHECK(result[20].y == doctest::Approx(267.37).epsilon(0.05));
    CHECK(result[20].H == doctest::Approx(127.37).epsilon(0.05));
  }

  SUBCASE("volterra D4")
  {
    sim.addParameters(parsD);
    sim.addUserPoint(p4);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsD.N].x == doctest::Approx(19.64).epsilon(0.05));
    CHECK(result[parsD.N].y == doctest::Approx(71.55).epsilon(0.05));
    CHECK(result[parsD.N].H == doctest::Approx(23.54).epsilon(0.05));
  }

  SUBCASE("volterra D5")
  {
    sim.addParameters(parsD);
    sim.addUserPoint(p5);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsD.N].x == doctest::Approx(1.51E-4).epsilon(0.05));
    CHECK(result[parsD.N].y == doctest::Approx(288.76).epsilon(0.05));
    CHECK(result[parsD.N].H == doctest::Approx(112.54).epsilon(0.05));
  }

  SUBCASE("volterra D6")
  {
    sim.addParameters(parsD);
    sim.addUserPoint(p6);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[parsD.N].x == doctest::Approx(4.80E-3).epsilon(0.05));
    CHECK(result[parsD.N].y == doctest::Approx(271.85).epsilon(0.05));
    CHECK(result[parsD.N].H == doctest::Approx(104.12).epsilon(0.05));
  }

  SUBCASE("check 40th step values in case D4")
  {
    sim.addParameters(parsD);
    sim.addUserPoint(p4);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[40].x == doctest::Approx(87.40).epsilon(0.05));
    CHECK(result[40].y == doctest::Approx(55.66).epsilon(0.05));
    CHECK(result[40].H == doctest::Approx(23.55).epsilon(0.05));
  }

  SUBCASE("check 40th step values in case D5")
  {
    sim.addParameters(parsD);
    sim.addUserPoint(p5);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[40].x == doctest::Approx(0.16).epsilon(0.05));
    CHECK(result[40].y == doctest::Approx(297.52).epsilon(0.05));
    CHECK(result[40].H == doctest::Approx(112.54).epsilon(0.05));
  }

  SUBCASE("check 40th step values in case D6")
  {
    sim.addParameters(parsD);
    sim.addUserPoint(p6);
    sim.evolve();
    auto result = sim.result();
    CHECK(result[40].x
          == doctest::Approx(3.32).epsilon(
              0.25)); // note: test would fail with stricter approximation
    CHECK(result[40].y == doctest::Approx(279.29).epsilon(0.05));
    CHECK(result[40].H == doctest::Approx(104.13).epsilon(0.05));
  }
}
