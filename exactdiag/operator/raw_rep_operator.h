#pragma once
#include "../global.h"

#include "../hilbertspace.h"

template <typename ...QNS>
class System;

template <typename _Scalar, typename ... QNS>
class RawRepOperator;

//! @class RawRepOperator
//! @brief RawRepOperator is the internal implementation of an operator. 
//!
//! @tparam _Scalar Scalar type
//! @tparam QNS List of quantum numbers (U(1))
template <typename _Scalar, typename ... QNS>
class RawRepOperator : public GenericOperator<_Scalar>
{
 public:
  using Scalar = _Scalar;
  using SystemType = System<QNS...>;

  //! Constructor
  //! @tparam RepSize Length of the binary representation
  //! @tparam SiteSize Number of sites (containts fermion parity, for example.)
  template <size_t RepSize, size_t SiteSize>
  RawRepOperator(const SystemType& sys,
                 PureOperator<Scalar, RepSize, SiteSize> const & po)
  : system_(sys)
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
        //size_t mdig = sys.site(idx_site).n_digit();
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


  //!< Getters
  std::vector<size_t> const & mask()     const { return mask_; }
  std::vector<size_t> const & row()      const { return row_; }
  std::vector<size_t> const & col()      const { return col_; }
  std::vector<size_t> const & fp_mask()  const { return fp_mask_; }
  std::vector<size_t> const & fp_row()   const { return fp_row_; }
  std::vector<size_t> const & fp_col()   const { return fp_col_; }
  std::vector<size_t> const & fp_check() const { return fp_check_; }
  Scalar coefficient() const { return coefficient_; }
  //!>


  //! conform
  //! @brief Check if the change of all quantum numbers is zero.
  bool conform(const System<QNS...>& system) const {
    auto row_quantum_number = std::make_tuple(QNS(0)...);
    auto col_quantum_number = std::make_tuple(QNS(0)...);
    size_t n_site = system.n_site();

    for (size_t i_site = 0 ; i_site < n_site ; ++i_site) {
      if (mask_[i_site]) {
        auto site = system.site(i_site);
        auto row_state = site.state(row_[i_site]);
        auto col_state = site.state(col_[i_site]);
        row_quantum_number = elementwise(row_quantum_number)
                             + elementwise(row_state.quantum_number());
        col_quantum_number = elementwise(col_quantum_number)
                             + elementwise(col_state.quantum_number());
      }
    }
    return row_quantum_number == col_quantum_number;
  }

 private:
  const SystemType& system_;
  std::vector<size_t> mask_, row_, col_;
  std::vector<size_t> fp_mask_, fp_row_, fp_col_, fp_check_;
  Scalar coefficient_;
};
