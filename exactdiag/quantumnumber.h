#pragma once

#include <cinttypes>
#include <tuple>
#include <type_traits>
#include <typeinfo>


namespace detail {
template<bool B, class T = void>
using enable_if_t = typename std::enable_if<B, T>::type;
}



template <typename _ValueType>
class QuantumNumber {
 public:
  using ValueType = _ValueType;

  static constexpr const char * name = "GenericQuantumNumber";

  QuantumNumber(const ValueType& value) : value_(value) { }
  const ValueType & value() const { return value_; }
 protected:
  ValueType value_;
};



template <typename QN,
          detail::enable_if_t<
                  std::is_base_of<QuantumNumber<typename QN::ValueType>, QN>::value,
                  int> = 0>
QN operator+(const QN& q1, const QN& q2)
{
  return QN(q1.value() + q2.value());
}

template <typename QN,
        detail::enable_if_t<
                std::is_base_of<QuantumNumber<typename QN::ValueType>, QN>::value,
                int> = 0>
QN operator-(const QN& q1, const QN& q2)
{
  return QN(q1.value() - q2.value());
}

template <typename QN,
        detail::enable_if_t<
                std::is_base_of<QuantumNumber<typename QN::ValueType>, QN>::value,
                int> = 0>
bool operator==(const QN& q1, const QN& q2)
{
  return (q1.value() == q2.value());
}


template <typename QN,
        detail::enable_if_t<
                std::is_base_of<QuantumNumber<typename QN::ValueType>, QN>::value,
                int> = 0>
bool operator<(const QN& q1, const QN& q2)
{
  return (q1.value() < q2.value());
}

template <typename QN,
        detail::enable_if_t<
                std::is_base_of<QuantumNumber<typename QN::ValueType>, QN>::value,
                int> = 0>
bool operator>(const QN& q1, const QN& q2)
{
  return (q1.value() > q2.value());
}

template <typename QN,
        detail::enable_if_t<
                std::is_base_of<QuantumNumber<typename QN::ValueType>, QN>::value,
                int> = 0>
bool operator<=(const QN& q1, const QN& q2)
{
  return (q1.value() <= q2.value());
}

template <typename QN,
        detail::enable_if_t<
                std::is_base_of<QuantumNumber<typename QN::ValueType>, QN>::value,
                int> = 0>
bool operator>=(const QN& q1, const QN& q2)
{
  return (q1.value() >= q2.value());
}

template <typename QN,
          detail::enable_if_t<
                  std::is_base_of<QuantumNumber<typename QN::ValueType>, QN>::value,
                  int> = 0>
std::ostream & operator<<(std::ostream& os, const QN& q)
{
  return (os << QN::name << "(" << q.value() << ")");
}



class Charge : public QuantumNumber<std::int64_t>
{
 public:
  static constexpr const char * name = "Charge";
  Charge() : QuantumNumber<std::int64_t>(0) { }
  Charge(std::int64_t value) : QuantumNumber<std::int64_t>(value) { }
};

class Spin : public QuantumNumber<std::int64_t>
{
 public:
  static constexpr const char * name = "Spin";
  Spin() : QuantumNumber<std::int64_t>(0) { }
  Spin(std::int64_t value) : QuantumNumber<std::int64_t>(value) { }
};

