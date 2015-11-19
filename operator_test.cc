#include "operator.h"
#include "quantumnumber.h"
#include "hilbertspace.h"
#include "tuple_tools.h"

int main(int argc, char** argv)
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
    auto const & site = system.site(i_site);
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

  return 0;
}
