#pragma once
#include <cassert>
#include <cmath>
#include <bitset>
#include <tuple>
#include <iostream>
#include <limits>
#include <vector>

template <typename _Scalar>
class GenericOperator {
 public:
  using Scalar = _Scalar;
};

template <typename _Scalar, size_t _RepSize, size_t _SiteSize>
class PureOperator;

template <typename _Scalar>
class RawRepOperator;

#include "pure_operator.h"
#include "raw_rep_operator.h"