#pragma once
#include "../global.h"

#include "../hilbertspace.h"

template <typename ...QNS>
class System;


template <typename _Scalar>
class RawRepOperator;


template <typename _Scalar>
class RawRepOperator : public GenericOperator<_Scalar>
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
