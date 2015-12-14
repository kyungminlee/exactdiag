#pragma once
#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include "quantumnumber.h"
#include "tuple_stream.h"
#include "tuple_tools.h"

#include "operator.h"

#define DEBUGRUN(X)

template <typename ... QNS>
class State {
 public:
  using QuantumNumberTuple = std::tuple<QNS...>;

  State(bool fermion_parity, const QNS & ... quantum_number)
          : name_("")
          , fermion_parity_(fermion_parity)
          , quantum_number_(quantum_number...)
  {
  }

  State(const char* name, bool fermion_parity, const QNS & ... quantum_number)
          : name_(name)
          , fermion_parity_(fermion_parity)
          , quantum_number_(quantum_number...)
  {
  }

  void display(std::ostream& os = std::cout, std::string prefix = "") const
  {
    os << prefix << "State" << std::endl;
    if (name_.size() > 0) { os << prefix << "  name          : " << name_ << std::endl; }
    os << prefix << "  fermion_parity: " << fermion_parity_ << std::endl;
    os << prefix << "  quantum_number: " << quantum_number_ << std::endl;
  }

  bool fermion_parity() const { return fermion_parity_; }
  const QuantumNumberTuple& quantum_number() const { return quantum_number_; }
  const QuantumNumberTuple& QN() const { return quantum_number_; }

 private:
  std::string name_;
  bool fermion_parity_;
  QuantumNumberTuple quantum_number_;
};


template <typename ...QNS>
class Site {
 public:
  using StateType = State<QNS...>;
  using QuantumNumberTuple = std::tuple<QNS...>;

  template <typename ... Args>
  Site(const StateType& state, Args... rest)
  {
    add_state(state, rest...);
  };

  Site& add_state() { return *this; }

  template <typename ...Args>
  Site& add_state(const StateType& state, Args ... args)
  {
    states_.push_back(state);
    return add_state(args...);
  }

  void display(std::ostream& os = std::cout, std::string prefix = "") const {
    os << prefix << "Site" << std::endl;
    for (auto const & state : states_ ) {
      state.display(os, prefix + std::string("  "));
    }
  }

  size_t n_state() const {
    return states_.size();
  }

  size_t n_digit() const {
    return static_cast<size_t>(std::ceil(::log2(states_.size())));
  }

  QuantumNumberTuple max_quantum_number() const {
    QuantumNumberTuple ret = states_[0].quantum_number();
    for (size_t i = 1 ; i < states_.size() ; ++i) {
      ret = elmax(ret, states_[i].quantum_number());
    }
    return ret;
  }

  QuantumNumberTuple min_quantum_number() const {
    QuantumNumberTuple ret = states_[0].quantum_number();
    for (size_t i = 1 ; i < states_.size() ; ++i) {
      ret = elmin(ret, states_[i].quantum_number());
    }
    return ret;
  }

  const StateType & state(size_t idx_state) const {
    assert(idx_state < states_.size());
    return states_[idx_state];
  }

 private:
  std::vector<StateType> states_;
};


template <size_t _RepSize, size_t _SiteSize, typename ... QNS>
class BasisIterator;

template <typename ...QNS>
class System {
 public:
  using StateType = State<QNS...>;
  using SiteType = Site<QNS...>;
  using QuantumNumberTuple = std::tuple<QNS...>;

  template<typename ...Args>
  System(Args ... sites) {
    add_site(sites...);
  }

  System & add_site() { return *this; }

  template<typename ... Args>
  System & add_site(const SiteType &site, Args ... args) {
    sites_.push_back(site);
    return add_site(args...);
  }


  SiteType & site(size_t idx_site) {
    assert(idx_site < n_site());
    return sites_[idx_site];
  }

  const SiteType & site(size_t idx_site) const {
    assert(idx_site < n_site());
    return sites_[idx_site];
  }

  void display(std::ostream &os = std::cout, std::string prefix = "") const {
    os << prefix << "System" << std::endl;
    for (auto const &site : sites_) {
      site.display(os, prefix + std::string("  "));
    }
  }

  size_t n_digit() const {
    size_t n_dig = 0;
    for (auto const &site : sites_) {
      n_dig += site.n_digit();
    }
    return n_dig;
  }

  size_t start_digit(size_t idx_site) const {
    assert(idx_site < n_site());
    size_t n_dig = 0;
    for (size_t j = 0; j < idx_site; ++j) {
      n_dig += sites_[j].n_digit();
    }
    return n_dig;
  }

  template<size_t RepSize>
  std::bitset<RepSize> mask_digit(size_t idx_site) const {
    assert(idx_site < n_site());
    size_t sdig = start_digit(idx_site);
    size_t mdig = sites_[idx_site].n_digit();
    std::bitset<RepSize> mask;
    for (size_t i = sdig; i < (sdig + mdig); ++i) {
      mask.set(i, true);
    }
    return mask;
  }

  template<size_t RepSize>
  void mask_digit(size_t idx_site, std::bitset<RepSize>& mask) const {
    assert(idx_site < n_site());
    size_t sdig = start_digit(idx_site);
    size_t mdig = sites_[idx_site].n_digit();
    mask.reset();
    for (size_t i = sdig; i < (sdig + mdig); ++i) {
      mask.set(i, true);
    }
  }


  // not inclusive
  QuantumNumberTuple max_quantum_number(size_t i_site) const {
    assert(i_site < sites_.size());
    QuantumNumberTuple c;
    for (size_t j = 0; j < i_site; ++j) {
      c = elementwise(c) + elementwise(sites_[j].max_quantum_number());
    }
    return c;
  }


  // not inclusive
  QuantumNumberTuple min_quantum_number(size_t i_site) const {
    assert(i_site < sites_.size());
    QuantumNumberTuple c;
    for (size_t j = 0; j < i_site; ++j) {
      c = elementwise(c) + elementwise(sites_[j].min_quantum_number());
    }
    return c;
  }


  QuantumNumberTuple max_quantum_number() const {
    QuantumNumberTuple c;
    for (size_t j = 0; j < n_site(); ++j) {
      c = elementwise(c) + elementwise(sites_[j].max_quantum_number());
    }
    return c;
  }


  QuantumNumberTuple min_quantum_number() const {
    QuantumNumberTuple c;
    for (size_t j = 0; j < n_site(); ++j) {
      c = elementwise(c) + elementwise(sites_[j].min_quantum_number());
    }
    return c;
  }


  template <typename PureOperatorType>
  bool comply(const PureOperatorType& op) {
    if (PureOperatorType::RepSize < n_digit()) { return false; }
    if (PureOperatorType::SiteSize < n_site()) { return false; }

  }


  template<typename Scalar, size_t RepSize, size_t SiteSize = RepSize>
  PureOperator<Scalar, RepSize, SiteSize>
  get_operator(size_t idx_site, size_t idx_rowstate, size_t idx_colstate) const {
    assert(n_digit() <= RepSize);
    assert(n_site() <= SiteSize);

    size_t nd = n_digit();
    size_t ns = sites_.size();

    assert(idx_site < ns);
    auto const &site = sites_[idx_site];
    assert(idx_rowstate < site.n_state());
    assert(idx_colstate < site.n_state());

    std::bitset<RepSize> r(idx_rowstate);
    std::bitset<RepSize> c(idx_colstate);

    r <<= start_digit(idx_site);
    c <<= start_digit(idx_site);
    std::bitset<RepSize> m = mask_digit<RepSize>(idx_site);

    auto row_parity = site.state(idx_rowstate).fermion_parity();
    auto col_parity = site.state(idx_colstate).fermion_parity();

    std::bitset<SiteSize> fm(1);
    fm <<= idx_site;

    std::bitset<SiteSize> fpr;
    if (row_parity) { fpr = fm; }
    std::bitset<SiteSize> fpc;
    if (col_parity) { fpc = fm; }
    std::bitset<SiteSize> fc;
    if (row_parity != col_parity) { // which means odd number of fermion operators.
      for (size_t j = 0; j < idx_site; ++j) {
        fc.set(j, true);
      }
    }

    Scalar v(1);
    return PureOperator<Scalar, RepSize, SiteSize>(m, r, c, fm, fpr, fpc, fc, v);
  }


  template<size_t RepSize, size_t SiteSize>
  std::tuple<std::bitset<RepSize>, std::bitset<SiteSize> >
  get_state_representation(size_t idx_site, size_t idx_state) const {

    assert(n_digit() <= RepSize);
    assert(n_site() <= SiteSize);

    size_t nd = n_digit();
    size_t ns = sites_.size();

    assert(idx_site < ns);
    auto const &site = sites_[idx_site];
    assert(idx_state < site.n_state());

    std::tuple<std::bitset<RepSize>, std::bitset<SiteSize> > ret;

    std::get<0>(ret) = idx_state;
    std::get<0>(ret) <<= start_digit(idx_site);

    if (site.state(idx_state).fermion_parity()) {
      std::get<1>(ret) = 1;
      std::get<1>(ret) <<= idx_site;
    }
    return ret;
  }

  template<size_t RepSize = 64, size_t SiteSize = RepSize>
  BasisIterator<RepSize, SiteSize, QNS...> cbegin(QNS... qns) const {
    return BasisIterator<RepSize, SiteSize, QNS...>(*this, qns...);
  };

  template<size_t RepSize = 64, size_t SiteSize = RepSize>
  BasisIterator<RepSize, SiteSize, QNS...> cend(QNS... qns) const {
    return BasisIterator<RepSize, SiteSize, QNS...>(BasisIterator<RepSize,SiteSize,QNS...>::InvalidIterator() );
  };



  size_t n_site() const { return sites_.size(); }

 private:
  std::vector<SiteType> sites_;
}; // class System


/*
 *    BI(S, Q) -> BI(S, Q) -> BI(S, Q) -> null
 *       0
 *
 *  for (alphabet A at S) {
 *    sub = get(S+1, Q(A))
 *    for (sub = begin(S+1, Q(A)) ; sub != end(S+1, Q(A)) ; ++sub)
 *    {
 *       ACT( [ A, *sub] )
 *    }
 *  }
 *
 */

template <size_t _RepSize, size_t _SiteSize, typename ... QNS>
class BasisIterator
{
 public:
  static const size_t RepSize = _RepSize;
  static const size_t SiteSize = _SiteSize;

  class InvalidIterator { };

  typedef System<QNS...> SystemType;
  using QuantumNumberTuple = typename SystemType::QuantumNumberTuple;

  using ValueType = std::tuple<std::bitset<RepSize>, std::bitset<SiteSize>>;
  using ConstReferenceType = const ValueType &;

  BasisIterator(InvalidIterator) : valid_(false) { }

  bool valid() const { return valid_; }

  BasisIterator(const SystemType& system,
                QNS ... quantum_number)
          : system_(system), quantum_number_(quantum_number...)
          , current_state_(system_.n_site())
          , cumulative_quantum_number_(system_.n_site()) // all zero at the beginning.
          , valid_(true)
  {
    assert(system_.n_site() > 0);
    cumulative_quantum_number_[system_.n_site()-1] = std::make_tuple(QNS(0)...);
    if (!set_first_rec(system_.n_site()-1)) {
      //throw std::logic_error("Dictionary empty");
      valid_ = false;
    }

  }

  template <typename T1, typename T2>
  using first_t = T1;

  static constexpr
  std::tuple<first_t<bool,QNS>...> all_true_quantum_number() {
    return std::make_tuple(first_t<bool,QNS>(true)...);
  };

  struct Tracker {
    Tracker(const char* in, const char* out)
    : in_(in), out_(out)
    {
      std::cout << in_ << std::endl;
    }
    ~Tracker(){
      std::cout << out_ << std::endl;
    }
    std::string in_, out_;
  };

  bool set_first_rec(size_t level) {
    assert(valid_);
    /*
    std::stringstream inmsg, outmsg;
    inmsg << "BEGIN set_first_rec(" << level << ")";
    outmsg << "END set_first_rec(" << level << ")";
    Tracker tt(inmsg.str().c_str(), outmsg.str().c_str());
    */
    DEBUGRUN(std::cout << "set_first_rec(" << level << ") : ";)
    DEBUGRUN(for(auto c : current_state_) { std::cout << c << ", ";} std::cout << std::endl;)

    if (level == 0) {
      auto const & site = system_.site(level);
      for (current_state_[level] = 0; current_state_[level] < site.n_state() ; ++current_state_[level]) {
        DEBUGRUN(std::cout << "LOOP(" << level << ") : " << current_state_[level] << std::endl;)
        auto const & state = site.state(current_state_[level]);
        auto qn = elementwise(cumulative_quantum_number_[level])
                  + elementwise(state.quantum_number());
        if (qn == quantum_number_) { return true; }
      }
      return false;
    } else if (level >= system_.n_site()) {
      throw std::logic_error("level too high");
    } else { // 0 < level < n_site
      auto const & site = system_.site(level);
      for (current_state_[level] = 0; current_state_[level] < site.n_state() ; ++current_state_[level]) {
        DEBUGRUN(std::cout << "LOOP(" << level << ") : " << current_state_[level] << std::endl;)
        auto const &state = site.state(current_state_[level]);
        auto qn = elementwise(cumulative_quantum_number_[level])
                  + elementwise(state.quantum_number());
        {
          auto min_qn = elementwise(qn)
                        + elementwise(system_.min_quantum_number(level));
          auto test_min = (elementwise(min_qn) <= elementwise(quantum_number_));
          DEBUGRUN(std::cout << "Testing Min QN: " << min_qn << " vs. " << quantum_number_ << " => "<< (all(test_min))<< std::endl;)
          if (!all(test_min)) { continue; }
        }
        {
          auto max_qn = elementwise(qn)
                        + elementwise(system_.max_quantum_number(level));
          auto test_max = (elementwise(quantum_number_) <= elementwise(max_qn));
          DEBUGRUN(std::cout << "Testing Max QN: " << max_qn << " vs. " << quantum_number_ << " => "<< (all(test_max))<< std::endl;)
          if (!all(test_max)) { continue; }
        }
        cumulative_quantum_number_[level-1] = qn;
        if (set_first_rec(level-1)) { return true; }
      }
      return false;
    }
  }


  BasisIterator& operator++()
  {
    if (!valid_) { return *this; }

    if( next(system_.n_site()-1)) {
      return *this;
    } else {
      valid_ = false;
      return *this;
    }
  }

  bool next(size_t level)
  {
    assert(valid_);
    DEBUGRUN(std::cout << "next(" << level << ") : ";)
    DEBUGRUN(for(auto c : current_state_) { std::cout << c << ", ";} std::cout << std::endl;)

    if (level == 0) {
      auto const & site = system_.site(level);
      for (++current_state_[level]; current_state_[level] < site.n_state() ; ++current_state_[level]) {
        auto const & state = site.state(current_state_[level]);
        auto qn = elementwise(cumulative_quantum_number_[level])
                  + elementwise(state.quantum_number());
        if (qn == quantum_number_) { return true; }
      }
      return false;
    } else if (level >= system_.n_site()) {
      throw std::logic_error("level too high");
    } else { // 0 < level < n_site
      if (next(level-1)) { return true; }
      auto const & site = system_.site(level);
      for (++current_state_[level]; current_state_[level] < site.n_state() ; ++current_state_[level]) {
        auto const &state = site.state(current_state_[level]);

        auto qn = elementwise(cumulative_quantum_number_[level])
                  + elementwise(state.quantum_number());
        {
          auto min_qn = elementwise(qn)
                        + elementwise(system_.min_quantum_number(level));
          auto test_min = (elementwise(min_qn) <= elementwise(quantum_number_));
          if (!all(test_min)) { continue; }
        }
        {
          auto max_qn = elementwise(qn)
                        + elementwise(system_.max_quantum_number(level));
          auto test_max = (elementwise(quantum_number_) <= elementwise(max_qn));
          if (!all(test_max)) { continue; }
        }
        cumulative_quantum_number_[level-1] = qn;
        if (set_first_rec(level-1)) { return true; }
      }
      return false;
    }
  }

  bool operator==(const BasisIterator& iter) const {
    if (! (valid_ || iter.valid_) ) { return true; } // both non-valid.

    return false; // TODO
  }

  const std::vector<size_t>& operator*() const {
    //if( !valid_) { throw std::runtime_error("trying to access invalid iterator"); }
    return current_state_;
  }

  std::tuple<std::bitset<RepSize>, std::bitset<SiteSize>>
  get() const
  {
    std::tuple<std::bitset<RepSize>, std::bitset<SiteSize> > ret;
    for (size_t i = 0 ; i < system_.n_site() ; ++i) {
      std::bitset<RepSize> vr(current_state_[i]);
      std::bitset<RepSize> vf(system_.site(i).state(current_state_[i]).fermion_parity());
      vr <<= system_.start_digit(i);
      vf <<= i;
      //auto func = system_.get_state_representation<RepSize, SiteSize>;
      //auto v = func(i, current_state_[i]);
      //auto v = system_.get_state_representation<RepSize, SiteSize>(i, current_state_[i]);
      std::get<0>(ret) |= vr;
      std::get<1>(ret) |= vf;
    }
    return ret;
  };

 private:
  const SystemType& system_;
  const QuantumNumberTuple quantum_number_;

  std::vector<QuantumNumberTuple> cumulative_quantum_number_;
  std::vector< size_t > current_state_;
  bool valid_;
  //std::tuple< std::bitset<RepSize>, std::bitset<SiteSize> > current_;
};
