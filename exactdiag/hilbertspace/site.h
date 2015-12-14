#pragma once
#include "../global.h"

#include "../utility/tuple_stream.h"
#include "../utility/tuple_tools.h"

#include "state.h"

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

