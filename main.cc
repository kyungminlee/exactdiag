#include <unordered_map>

#include "operator.h"
#include "hilbertspace.h"
#include "utility/tuple_tools.h"


int main_spinfermion(int argc, char** argv)
{
  static const size_t RepSize = 9;
  static const size_t SiteSize = 9;

  using namespace std;
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
  //system.display();


  MixedOperator<double, RepSize, SiteSize> ann, cre, hop;
  for (size_t i_site = 0 ; i_site < system.n_site() ; ++i_site) {
    //auto const &site = system.site(i_site);
    auto c = system.get_operator<double, RepSize, SiteSize>(i_site, 0, 1);
    auto cdag = system.get_operator<double, RepSize, SiteSize>(i_site, 1, 0);
    ann.add(c);
    cre.add(cdag);
    RawRepOperator<double, Charge, Spin> raw_ann(system, c);
    c.display();
    raw_ann.display();
    cout << endl;
  }

  //cout << "===== Creation =====" << endl;
  //cre.display();
  //cout << "===== Annihilation =====" << endl;
  //ann.display();

  for (size_t i_cell = 0 ; i_cell < 3 ; ++i_cell) {
    for (size_t i_site = 3*i_cell + 1; i_site < 3*i_cell + 3 ; ++i_site) {
      size_t j_site = (i_site + 3) % system.n_site();
      hop.add(-0.5
              * system.get_operator<double, RepSize, SiteSize>(i_site, 1, 0)
              * system.get_operator<double, RepSize, SiteSize>(j_site, 0, 1));
      hop.add(-0.5
              * system.get_operator<double, RepSize, SiteSize>(j_site, 1, 0)
              * system.get_operator<double, RepSize, SiteSize>(i_site, 0, 1));
    }
  }
  //cout << "===== Hoppings =====" << endl;
  //hop.display();
  //return 0;

  //cout << "===== Basis Generation =====" << endl;

  vector<tuple<bitset<RepSize>, bitset<SiteSize>>> basis_list;
  unordered_map<bitset<RepSize>, size_t> basis_map;
  for (auto iter = system.cbegin<RepSize, SiteSize>(Charge(3), Spin(0)) ; iter.valid() ; ++iter) {
    basis_list.push_back(iter.get());
    //cout << "BASIS GENERATED : " << iter.get() << endl;
  }
  for (size_t i = 0 ; i < basis_list.size() ; ++i) {
    basis_map[std::get<0>(basis_list[i])] = i;
    cout << i << "\t" << std::get<0>(basis_list[i]) << endl;
  }


  for (size_t i_basis = 0 ; i_basis < basis_list.size() ; ++i_basis) {
    auto const & bvec_fvec = basis_list[i_basis];
    auto row = hop.apply(std::get<0>(bvec_fvec), std::get<1>(bvec_fvec));
    for (auto const & r : row) {
      auto match_iter = basis_map.find(std::get<0>(r));
      //cout << std::get<0>(r) << "\t" << std::get<0>(bvec_fvec) << "\t" << std::get<2>(r) << endl;
      if (match_iter == basis_map.end()) {
        cout << "ERROR!" << endl;
      } else {
        cout << match_iter->second << "\t" << i_basis << "\t" << std::get<2>(r) << endl;
      }
    }
  }
  return 0;
}


#if 0
struct Generator
{
  template <size_t S1, size_t S2>
  std::tuple<std::bitset<S1>, std::bitset<S2>> myfun() const {
    std::tuple<std::bitset<S1>, std::bitset<S2>> ret;
    return ret;
  }

};

template<typename GeneratorType, size_t S1, size_t S2>
struct Wrapper
{
  std::tuple<std::bitset<S1>, std::bitset<S2>> myfun() const {
    return g.myfun<S1, S2>();
  }
  GeneratorType g;
};
#endif

int main(int argc, char** argv)
{

  //Wrapper<Generator, 3, 4> wrap;
  //std::cout << wrap.myfun() << std::endl;
  return main_spinfermion(argc, argv);
  return 0;
}

int main_test(int argc, char** argv)
{
  using namespace std;

  cout << "===== Operator Test =====" << endl;
  PureOperator<double, 4, 4> op;
  op.display();

  cout << "===== Quantum Number Test =====" << endl;
  Charge q1(1);
  Charge q2(2);


  cout << q1 << endl;
  cout << q2 << endl;
  cout << (q1+q2) << endl;
  cout << (q1-q2) << endl;

  cout << "===== State Test =====" << endl;

  State<Charge,Spin> emp(false,  Charge(0), Spin( 0));
  State<Charge,Spin> occ(false,  Charge(1), Spin( 0));
  State<Charge,Spin> up( true,   Charge(1), Spin( 1));
  State<Charge,Spin> dn( true,   Charge(1), Spin(-1));
  State<Charge,Spin> updn(false, Charge(2), Spin( 0));

  emp.display();
  up.display();
  dn.display();
  updn.display();

  cout << "===== Site Test =====" << endl;
  Site<Charge, Spin> site(emp, occ);//, dn, updn);
  site.display();
  cout << "max QN  : " << site.max_quantum_number() << endl;
  cout << "min QN  : " << site.min_quantum_number() << endl;
  cout << "n_digit : " << site.n_digit() << endl;

  cout << "===== System Test =====" << endl;
  System<Charge, Spin> system(site, site, site, site, site, site);
  system.display();
  cout << "max QN  : " << system.max_quantum_number() << endl;
  cout << "min QN  : " << system.min_quantum_number() << endl;
  cout << "n_digit : " << system.n_digit() << endl;

  cout << "===== System Operator Test =====" << endl;
  MixedOperator<double, 6, 6> ann, cre;
  for (size_t i_site = 0 ; i_site < system.n_site() ; ++i_site) {
    //auto const & site = system.site(i_site);
    auto c    = system.get_operator<double, 6, 6>(i_site, 0, 1);
    auto cdag = system.get_operator<double, 6, 6>(i_site, 1, 0);
    ann.add(c);
    cre.add(cdag);
  }

  cout << "===== Annihilation =====" << endl;
  ann.display();
  cout << "===== Creation =====" << endl;
  cre.display();

  cout << "===== Hopping =====" << endl;
  (cre.term(4) * ann.term(1)).display();

  cout << "===== Basis Generation =====" << endl;

  for (auto iter = system.cbegin(Charge(3), Spin(0)) ; iter.valid() ; ++iter) {

    for (auto f : *iter) { cout << f << "\t"; } cout << endl;

  }
  return 0;
}
