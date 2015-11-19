#include <tuple>
#include <iostream>
#include <functional>
#include <utility>
#include "tuple_stream.h"
#include "tuple_tools.h"

int main()
{
  using namespace std;
  const tuple<int, int, int> a(100, 200,20);
  tuple<int, int, int> b(100,1000,10);
  auto op = make_tuple(plus<int>(), plus<int>(), plus<int>());
  cout << "a = " << a << endl;
  cout << "b = " << b << endl;
  auto c = elementwise(a).apply_binary(op, elementwise(b));

#if 0
  cout << "a +  b = " << (elementwise(a) +  elementwise(b)) << endl;
  cout << "a -  b = " << (elementwise(a) -  elementwise(b)) << endl;
  cout << "a *  b = " << (elementwise(a) *  elementwise(b)) << endl;
  cout << "a /  b = " << (elementwise(a) /  elementwise(b)) << endl;
  cout << "a == b = " << (elementwise(a) == elementwise(b)) << endl;
  cout << "a != b = " << (elementwise(a) != elementwise(b)) << endl;
  cout << "a <= b = " << (elementwise(a) <= elementwise(b)) << endl;
  cout << "a >= b = " << (elementwise(a) >= elementwise(b)) << endl;
  cout << "a <  b = " << (elementwise(a) <  elementwise(b)) << endl;
  cout << "a >  b = " << (elementwise(a) >  elementwise(b)) << endl;
  cout << "max(a,b)= " << elmax(elementwise(a), elementwise(b)) << endl;
  cout << "min(a,b)= " << elmin(elementwise(a), elementwise(b)) << endl;
  auto c = elmin(elementwise(a), elementwise(b));
  //auto c = (elementwise(a) > elementwise(b));
#endif
  cout << "a = " << a << endl;
  cout << "b = " << b << endl;
  //cout << "c = " << c << endl;
  //cout << "c.any() : " << any(c) << endl;
  //cout << "c.all() : " << all(c) << endl;
  return 0;
}

