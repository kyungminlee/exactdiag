#pragma once
#include "../global.h"

#include "../utility/tuple_stream.h"
#include "../utility/tuple_tools.h"

#include "state.h"

//! @class Site
//! @ingroup Hilbertspace_Module
//!
//! @brief Site defined by list of states.
//!
//! The %Site class represents an abstract site which can have a number of
//! different states.
//!
//! @tparam QNS List of U(1) quantum numbers
template <typename ...QNS>
class Site {
 public:
  using StateType = State<QNS...>;
  using QuantumNumberTuple = std::tuple<QNS...>;

  //! Constructor.
  //!
  //! Start with at least one state.
  //! @param state State
  template <typename ... Args>
  Site(const StateType& state, Args... rest)
  {
    add_state(state, rest...);
  };

private:
  Site& add_state() { return *this; }

public:

  //! Add possible states to the Site
  template <typename ...Args>
  Site& add_state(const StateType& state, Args ... args)
  {
    states_.push_back(state);
    return add_state(args...);
  }

  //! Prettyprint a %Site.
  void display(std::ostream& os = std::cout, std::string prefix = "") const {
    os << prefix << "Site" << std::endl;
    for (auto const & state : states_ ) {
      state.display(os, prefix + std::string("  "));
    }
  }

  //! Number of possible states of a site.
  size_t n_state() const {
    return states_.size();
  }

  //! Number of binary digits required to represent this %Site.
  size_t n_digit() const {
    return static_cast<size_t>(std::ceil(::log2(states_.size())));
  }

  //! Maximum possible U(1) quantum number of this %Site.
  //! @return Tuple of integers
  QuantumNumberTuple max_quantum_number() const {
    QuantumNumberTuple ret = states_[0].quantum_number();
    for (size_t i = 1 ; i < states_.size() ; ++i) {
      ret = elmax(ret, states_[i].quantum_number());
    }
    return ret;
  }

  //! Minimum possible U(1) quantum number of this %Site.
  //! @return Tuple of integers
  QuantumNumberTuple min_quantum_number() const {
    QuantumNumberTuple ret = states_[0].quantum_number();
    for (size_t i = 1 ; i < states_.size() ; ++i) {
      ret = elmin(ret, states_[i].quantum_number());
    }
    return ret;
  }

  //! Get a state with the given index.
  //! @return State
  const StateType & state(size_t idx_state) const {
    assert(idx_state < states_.size());
    return states_[idx_state];
  }

 private:
  std::vector<StateType> states_;
};
