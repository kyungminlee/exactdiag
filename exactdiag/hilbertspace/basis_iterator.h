#pragma once
#include "../global.h"

#include "system.h"

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

  BasisIterator(const SystemType& system,
                QNS ... quantum_number)
      : system_(system), quantum_number_(quantum_number...)
      , cumulative_quantum_number_(system_.n_site()) // all zero at the beginning.
      , current_state_(system_.n_site())
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

  bool valid() const { return valid_; }

 private:
  const SystemType& system_;
  const QuantumNumberTuple quantum_number_;

  std::vector<QuantumNumberTuple> cumulative_quantum_number_;
  std::vector< size_t > current_state_;
  bool valid_;
  //std::tuple< std::bitset<RepSize>, std::bitset<SiteSize> > current_;
};

