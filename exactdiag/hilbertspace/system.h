#pragma once
#include "../global.h"

#include "quantumnumber.h"
#include "state.h"
#include "site.h"
#include "operator.h"

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
    size_t nd = n_digit();
    size_t ns = sites_.size();

    assert(nd <= RepSize);
    assert(ns <= SiteSize);
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

