#pragma once
#include "../global.h"

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

