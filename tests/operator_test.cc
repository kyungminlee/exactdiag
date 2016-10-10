#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include <iostream>

#include "operator.h"
#include "hilbertspace.h"

TEST_CASE("Site test", "[site]") {
  State<Spin> su("SpinUp", false, Spin(1));
  State<Spin> s0("SpinZr", false, Spin(0));
  State<Spin> sd("SpinDn", false, Spin(-1));
  Site<Spin> spin_site1(su, s0, sd);

  //INFO( "Details of the created Site (a S = 1 spin site)" );
  spin_site1.display(Catch::cout());

  SECTION( "all member methods function properly") {
    REQUIRE(spin_site1.n_state() == 3);
    REQUIRE(spin_site1.n_digit() == 2);
    auto maxqn = spin_site1.max_quantum_number();
    auto minqn = spin_site1.min_quantum_number();
    REQUIRE(maxqn == std::make_tuple(Spin(1)));
    REQUIRE(minqn == std::make_tuple(Spin(-1)));
    REQUIRE(spin_site1.state(0) == su);
    REQUIRE(spin_site1.state(1) == s0);
    REQUIRE(spin_site1.state(2) == sd);
  }

  SECTION("add_state works properly") {
    Site<Spin> spin_site2(su);
    spin_site2.add_state(s0)
        .add_state(sd);
    REQUIRE(spin_site1 == spin_site2);

    Site<Spin> spin_site3(sd);
    spin_site3.add_state(s0).add_state(su);

    REQUIRE(spin_site1 != spin_site3);
  }
}

TEST_CASE("Ising model test", "[ising]") {
  //static const size_t RepSize = 16;
  //static const size_t SiteSize = 16;
  System<Spin> ising_system;
  {
    State<Spin> su("SpinUp", false, Spin(1));
    State<Spin> sd("SpinDn", false, Spin(-1));
    Site<Spin> spin_site(su, sd);
    for (size_t i = 0; i < 8; ++i) {
      ising_system.add_site(spin_site);
    }
  }

  {
    size_t count = 0;

    for (auto const & v : ising_system.subspace(Spin(1))) {
      ++count;
    }
    REQUIRE(count == 0);
  }

  {
    size_t count = 0;
    for (auto const & v : ising_system.subspace(Spin(2))) {
      size_t sum = 0;
      for (auto const & i : v) {
        sum += i;
      }
      REQUIRE(sum == 3);
      ++count;
    }
    REQUIRE(count == 56);
  }
}

TEST_CASE("Spin-fermion model test", "[spin-fermion]") {
  static const size_t RepSize = 9;
  static const size_t SiteSize = 9;
  System<Charge, Spin> system;
  {
    State<Charge, Spin> f0("FEm", false, Charge(0), Spin(0));
    State<Charge, Spin> fu("FUp", true, Charge(1), Spin(1));
    State<Charge, Spin> fd("FDn", true, Charge(1), Spin(-1));

    State<Charge, Spin> su("SUp", false, Charge(0), Spin(1));
    State<Charge, Spin> sd("SDn", false, Charge(0), Spin(-1));

    Site<Charge, Spin> spin_site(su, sd);
    Site<Charge, Spin> fup_site(f0, fu);
    Site<Charge, Spin> fdn_site(f0, fd);


    system.add_site(spin_site);
    system.add_site(fup_site);
    system.add_site(fdn_site);
    system.add_site(spin_site);
    system.add_site(fup_site);
    system.add_site(fdn_site);
    system.add_site(spin_site);
    system.add_site(fup_site);
    system.add_site(fdn_site);
  }

  MixedOperator<double, RepSize, SiteSize> ann, cre, hop;
  for (size_t i_site = 0 ; i_site < system.n_site() ; ++i_site) {
    auto const &site = system.site(i_site);
    auto c = system.get_operator<double, RepSize, SiteSize>(i_site, 0, 1);
    auto cdag = system.get_operator<double, RepSize, SiteSize>(i_site, 1, 0);
    ann.add(c);
    cre.add(cdag);
  }

  for (size_t i_cell = 0 ; i_cell < 3 ; ++i_cell) {
    for (size_t i_site = 3*i_cell + 1; i_site < 3*i_cell + 3 ; ++i_site) {
      size_t j_site = (i_site + 3) % system.n_site();
      hop.add( -0.5
               * system.get_operator<double, RepSize, SiteSize>(i_site, 1, 0)
               * system.get_operator<double, RepSize, SiteSize>(j_site, 0, 1));
      hop.add( -0.5
               * system.get_operator<double, RepSize, SiteSize>(j_site, 1, 0)
               * system.get_operator<double, RepSize, SiteSize>(i_site, 0, 1));
    }
  }


}