#pragma once
//!
//!

#include "../global.h"

//! @class SectorGenerator.
//!
//! @brief Generator of Sector defined by QuantumNumbers.
//!
//! @tparam _RepSize  Number of bits of representation (for basis representation).
//! @tparam _SiteSize  Number of sites (for site representation)
//! @tparam QuantumNumbers List of U(1) quantum numbers.
//!
//! A Sector class contains a vector of tuple of BitRep and BitSite, and unordered map of BitRep to ...
//! (forward and backward map between representation and index)
template <size_t _RepSize, size_t _SiteSize, typename ...QuantumNumbers>
class SectorGenerator {
 public:
  using SystemType = System<QuantumNumbers...>;
  static const size_t RepSize = _RepSize;
  static const size_t SiteSize = _SiteSize;
  using BitRep = std::bitset<RepSize>;
  using BitSite = std::bitset<SiteSize>;

 public:
  struct Sector
  {
    std::vector<std::tuple<BitRep, BitSite>> basis;
    std::unordered_map<BitRep, size_t> basismap;
  };

  //! Constructor with the given System.
  //! @param system %System
  SectorGenerator(const SystemType& system)
      : system_(system)
  {
  }

  //! @brief Generate a sector of the Hilbertspace with the given quantum numbers.
  //! @param qns List of quantum numbers.
  Sector generate(QuantumNumbers... qns) const
  {
    Sector sector;
    for (auto iter = system_. template cbegin<RepSize, SiteSize>(qns...);
         iter.valid();
         ++iter) {
      sector.basis.push_back(iter.get());
    }
    for (size_t i = 0; i < sector.basis.size(); ++i) {
      sector.basismap[std::get<0>(sector.basis[i])] = i;
    }
    return sector;
  }

 private:
  const SystemType& system_;
};
