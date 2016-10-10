#pragma once
#include "../global.h"

//! State
//! defined by fermion parity and quantum number.
//!
//! @tparam QNS List of U(1) quantum numbers.
template <typename ... QNS>
class State {
 public:
  using QuantumNumberTuple = std::tuple<QNS...>;

  //! Constructor
  //! @brief Construct a %State of a %Site with the given fermion parity and quantum numbers
  //! @param fermion_parity
  //! @param quantum_number
  State(bool fermion_parity, const QNS & ... quantum_number)
      : name_("")
      , fermion_parity_(fermion_parity)
      , quantum_number_(quantum_number...)
  {
  }

  //! Constructor
  //! @brief Construct a %State of a %Site with the given name, fermion parity and quantum numbers
  //! @param name
  //! @param fermion_parity
  //! @param quantum_number
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

