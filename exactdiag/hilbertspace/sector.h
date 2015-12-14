#pragma once
#include "../global.h"


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

  SectorGenerator(const SystemType& system)
      : system_(system)
  {
  }

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

