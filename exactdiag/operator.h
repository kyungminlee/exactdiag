#pragma once
#include <cassert>
#include <cmath>
#include <bitset>
#include <tuple>
#include <iostream>
#include <limits>
#include <vector>

class GenericOperator{ };
template <typename _Scalar, size_t _RepSize, size_t _SiteSize>
class PureOperator;

template <typename ...QNS>
class System;

template <typename _Scalar>
class RawRepOperator
{
 public:
  using Scalar = _Scalar;

  template <size_t RepSize, size_t SiteSize, typename ... QNS>
  RawRepOperator(const System<QNS...>& sys, PureOperator<Scalar, RepSize, SiteSize> const & po)
  {
    size_t ns = sys.n_site();
    size_t nd = sys.n_digit();
    assert(nd >= RepSize);
    assert(ns >= SiteSize);

    auto m = po.mask();
    for (size_t idx_site = 0 ; idx_site < ns ; ++idx_site) {
      std::bitset<RepSize> site_mask;
      sys.mask_digit(idx_site, site_mask);
      if ( (m & site_mask).any()) {
        assert((m & site_mask).count() == site_mask.count());
        size_t sdig = sys.start_digit(idx_site);
        size_t mdig = sys.site(idx_site).n_digit();
        auto r = (po.row() & site_mask) >> sdig;
        auto c = (po.col() & site_mask) >> sdig;
        mask_.push_back(idx_site);
        row_.push_back(r.to_ulong());
        col_.push_back(c.to_ulong());
      }

      std::bitset<SiteSize> site_fmask = 1;
      site_fmask <<= idx_site;
      if ((po.fp_mask() & site_fmask).any()) {
        auto r = (po.fp_row() & site_fmask) >> idx_site;
        auto c = (po.fp_col() & site_fmask) >> idx_site;
        fp_mask_.push_back(idx_site);
        fp_row_.push_back(r.to_ulong());
        fp_col_.push_back(c.to_ulong());
      }

      if ((po.fp_check() & site_fmask).any()) {
        fp_check_.push_back(idx_site);
      }
    }
    coefficient_ = po.coefficient();
  }

  void display(std::ostream& os = std::cout, std::string prefix = "") const {
    auto show = [&os, &prefix](const char* str, const std::vector<size_t>& vlist) -> void {
      os << prefix << str;
      bool first = true;
      for (auto v : vlist) { if (!first) { os << ", "; } os << v; first = false; }
      os << std::endl;
    };
    os << prefix << "RawRepOperator(" << this << ")" << std::endl;
    show(           "  mask        : ", mask_);
    show(           "  row         : ", row_);
    show(           "  col         : ", col_);
    show(           "  fp_mask     : ", fp_mask_);
    show(           "  fp_row      : ", fp_row_);
    show(           "  fp_col      : ", fp_col_);
    show(           "  fp_check    : ", fp_check_);
    os << prefix << "  coefficient : " << coefficient_ << std::endl;
  }

 private:
  std::vector<size_t> mask_, row_, col_;
  std::vector<size_t> fp_mask_, fp_row_, fp_col_, fp_check_;
  Scalar coefficient_;
};

template <typename _Scalar, size_t _RepSize, size_t _SiteSize>
class PureOperator : public GenericOperator
{
 public:
  static const size_t RepSize = _RepSize;
  static const size_t SiteSize = _SiteSize;

  using Scalar = _Scalar;
  using Rep = std::bitset<RepSize>;
  using SiteRep = std::bitset<SiteSize>;

  static_assert(SiteSize <= RepSize, "SiteSize cannot be greater than RepSize");

  PureOperator()
      : coefficient_( Scalar(0) )
  {
  }

  PureOperator(const std::bitset<RepSize> & mask,
               const std::bitset<RepSize> & row,
               const std::bitset<RepSize> & col,
               const std::bitset<SiteSize> & fp_mask,
               const std::bitset<SiteSize> & fp_row,
               const std::bitset<SiteSize> & fp_col,
               const std::bitset<SiteSize> & fp_check,
               const Scalar& coeff)
      : mask_(mask), row_(row), col_(col)
      , fp_mask_(fp_mask), fp_row_(fp_row), fp_col_(fp_col), fp_check_(fp_check)
      , coefficient_(coeff)
  {
    if ( (~mask & row).any() || (~mask & col).any() ) {
      throw std::domain_error("PureOperator(): Wrong overlap (m vs. r/c)");
    } else if ( (~fp_mask & fp_row).any() || (~fp_mask & fp_col).any()) {
      throw std::domain_error("PureOperator(): Wrong overlap (fm vs. fpr/fpc)");
    } else if ( (fp_check & fp_row).any() || (fp_check & fp_col).any()) {
      throw std::domain_error("PureOperator(): Wrong overlap (fc vs. fpr/fpc)");
    }
  }

  std::bitset<RepSize>  const & mask() const { return mask_; }
  std::bitset<RepSize>  const & row() const { return row_; }
  std::bitset<RepSize>  const & col() const { return col_; }
  std::bitset<SiteSize> const & fp_mask() const { return fp_mask_; }
  std::bitset<SiteSize> const & fp_row() const { return fp_row_; }
  std::bitset<SiteSize> const & fp_col() const { return fp_col_; }
  std::bitset<SiteSize> const & fp_check() const { return fp_check_; }
  Scalar const & coefficient() const { return coefficient_; }


  PureOperator operator*(Scalar v) const {
    if (std::abs(v) < std::numeric_limits<Scalar>::epsilon()) {
      return PureOperator();
    } else {
      return PureOperator(mask_, row_, col_, fp_mask_, fp_row_, fp_col_, fp_check_, coefficient_ * v);
    }
  }

  PureOperator operator*(const PureOperator& rhs) const {
    auto const & x = *this;
    auto const & y = rhs;
    
    auto ma = x.mask_ & y.mask_;
    if ((ma & x.col_) != (ma & y.row_) ) {
      return PureOperator();
    } else {
      auto mo = x.mask_ | y.mask_;
      auto mx = x.mask_ ^ ma;
      auto my = y.mask_ ^ ma;
      auto fma = x.fp_mask_ & y.fp_mask_;
      auto fmo = x.fp_mask_ | y.fp_mask_;
      auto fmx = x.fp_mask_ ^ fma;
      auto fmy = y.fp_mask_ ^ fma;

      assert( (fma & x.fp_col_) == (fma & y.fp_row_) );
      auto m = mo;
      auto r = x.row_ | (y.row_ & my);
      auto c = y.col_ | (x.col_ & mx);
      auto fm = fmo;
      auto fpr = x.fp_row_ | (y.fp_row_ & fmy);
      auto fpc = y.fp_col_ | (x.fp_col_ & fmx);

      auto parity_sign = [](const std::bitset<SiteSize>& bs) -> int {
        return ( bs.count() % 2 == 0 ) ? 1 : -1;
      };
      auto turn_off = [](std::bitset<SiteSize>& tgt, const std::bitset<SiteSize>& mask) -> void {
        tgt ^= (tgt & mask);
      };
      
      int sgn = 1;
      auto fc = x.fp_check_;
      sgn *= parity_sign(fc & y.fp_row_);
      turn_off(fc, y.fp_mask_);
      fc ^= y.fp_check_;

      sgn *= parity_sign(fc & fpc);
      turn_off(fc, fm);

      auto v = x.coefficient_ * y.coefficient_ * sgn;
      return PureOperator(m, r, c, fm, fpr, fpc, fc, v);
    }
  }

  bool match(const std::bitset<RepSize>& bvec) const {
    return ((bvec & mask_) == col_);
  }

  std::tuple<Rep, SiteRep, Scalar> apply(const Rep& bvec, const SiteRep& fvec) const {
    assert(match(bvec));
    auto fret = (fvec & ~fp_mask_) | (fp_row_);
    auto bret = (bvec & ~mask_) | (row_);
    auto sgn = ((fvec & fp_check_).count() % 2 == 0) ? 1 : -1;
    auto coeff = coefficient_ * sgn;
    return std::make_tuple(bret, fret, coeff);
  }

  void display(std::ostream& os = std::cout, std::string prefix = "") const {
    os << prefix << "PureOperator(" << this << ")" << std::endl;
    os << prefix << "  mask        : " << mask_ << std::endl;
    os << prefix << "  row         : " << row_ << std::endl;
    os << prefix << "  col         : " << col_ << std::endl;
    os << prefix << "  fp_mask     : " << fp_mask_ << std::endl;
    os << prefix << "  fp_row      : " << fp_row_ << std::endl;
    os << prefix << "  fp_col      : " << fp_col_ << std::endl;
    os << prefix << "  fp_check    : " << fp_check_ << std::endl;
    os << prefix << "  coefficient : " << coefficient_ << std::endl;
  }

 private:
  std::bitset<RepSize> mask_, row_, col_;
  std::bitset<SiteSize> fp_mask_, fp_row_, fp_col_, fp_check_;
  Scalar coefficient_;
};

template <typename Scalar, size_t RepSize, size_t SiteSize> inline
PureOperator<Scalar, RepSize, SiteSize> operator*(const Scalar & v,
                                                  const PureOperator<Scalar, RepSize, SiteSize> & op)
{
  if (std::abs(v) < std::numeric_limits<Scalar>::epsilon()) {
    return PureOperator<Scalar, RepSize, SiteSize>();
  } else {
    return PureOperator<Scalar, RepSize, SiteSize>(
            op.mask(), op.row(), op.col(),
            op.fp_mask(), op.fp_row(), op.fp_col(), op.fp_check(),
            v * op.coefficient());
  }
};



//===============


template <typename _Scalar, size_t _RepSize, size_t _SiteSize = _RepSize>
class MixedOperator : public GenericOperator
{
 public:
  static const size_t RepSize = _RepSize;
  static const size_t SiteSize = _SiteSize;

  using Scalar = _Scalar;
  using Rep = std::bitset<RepSize>;
  using SiteRep = std::bitset<SiteSize>;
  using PureOperatorType = PureOperator<Scalar, RepSize, SiteSize>;

  MixedOperator()
  {
  }

  template <typename T>
  MixedOperator& add(T&& op) {
    terms_.push_back(std::forward<T>(op));
    return *this;
  }

  const PureOperatorType & term(size_t i_term) const {
    assert(i_term < terms_.size());
    return terms_[i_term];
  }


  void display(std::ostream& os = std::cout, std::string prefix = "") const
  {
    os << prefix << "MixedOperator" << std::endl;
    for (auto const & term : terms_) {
      term.display(os, prefix + "  ");
    }
  }

  std::vector<std::tuple<Rep, SiteRep, Scalar>>
  apply(const Rep& bvec, const SiteRep& fvec) const {
    std::vector<std::tuple<Rep, SiteRep, Scalar>> ret;
    for (auto const & term : terms_) {
      if (term.match(bvec)) {
        ret.push_back(term.apply(bvec, fvec));
      }
    }
    return ret;
  }

 private:
  std::vector<PureOperatorType> terms_;
};