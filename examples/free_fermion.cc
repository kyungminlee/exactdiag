//
// Created by kmlee on 12/13/2015.
//

#include <cinttypes>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <kore/array/array.h>
#include "hilbertspace.h"
#include "operator.h"
#include <Eigen/Eigen>
#include <Eigen/SparseCore>


template <typename _Scalar>
class ChargeSpinSystemBuilder
{
 public:
  using Scalar = _Scalar;
  static const size_t RepSize = 64;
  static const size_t SiteSize = 64;

  using SystemType = System<Charge, Spin>;
  using PureOperatorType = PureOperator<double, RepSize, SiteSize>;
  using MixedOperatorType = MixedOperator<double, RepSize, SiteSize>;
  using StateType = State<Charge, Spin>;
  using SiteType = Site<Charge, Spin>;

  enum { kSpinUp = 0, kSpinDn = 1, kNumSpin = 2 };

  ChargeSpinSystemBuilder(int64_t nx, int64_t ny) {

  }


 private:

};


int main(int argc, char** argv)
{
  using namespace std;
  static const size_t RepSize = 32;
  static const size_t SiteSize = 32;

  using SystemType = System<Charge, Spin>;
  using PureOperatorType = PureOperator<double, RepSize, SiteSize>;
  using MixedOperatorType = MixedOperator<double, RepSize, SiteSize>;
  using StateType = State<Charge, Spin>;
  using SiteType = Site<Charge, Spin>;

  enum { kSpinUp = 0, kSpinDn = 1, kNumSpin = 2 };
  size_t nx = 4;
  size_t ny = 4;

  SystemType system;

  {
    StateType f0("FEm", false, Charge(0), Spin(0));
    StateType fu("FUp", true, Charge(1), Spin(1));
    StateType fd("FDn", true, Charge(1), Spin(-1));

    SiteType fup_site(f0, fu);
    SiteType fdn_site(f0, fd);

    for (size_t ix = 0 ; ix < nx ; ++ix) {
      for (size_t iy = 0; iy < ny; ++iy) {
        system.add_site(fup_site);
        system.add_site(fdn_site);
      }
    }
  }

  std::vector<PureOperatorType> ann_data, cre_data;
  for (size_t ix = 0 ; ix < nx ; ++ix) {
    for (size_t iy = 0 ; iy < ny ; ++iy) {
      for (size_t i_spin = 0 ; i_spin < 2 ; ++i_spin) {
        size_t i_site = (ix * ny + iy) * 2 + i_spin;
        ann_data.push_back(system.get_operator<double, RepSize, SiteSize>(i_site, 0, 1));
        cre_data.push_back(system.get_operator<double, RepSize, SiteSize>(i_site, 1, 0));
      }
    }
  }

  kore::array::Array<PureOperatorType, 3> ann(ann_data.data(), nx, ny, kNumSpin);
  kore::array::Array<PureOperatorType, 3> cre(cre_data.data(), nx, ny, kNumSpin);

  double t0 = -1.5;
  double t1 = 1.0;

  std::vector<double> vs  = {-t0, -t1, -t1, -t1, -t1};
  std::vector<int>    dxs = {  0,   0,   0,   1,  -1};
  std::vector<int>    dys = {  0,   1,  -1,   0,   0};
  size_t n_hop = vs.size();

  MixedOperatorType hop;

  for (size_t ix = 0 ; ix < nx ; ++ix) {
    for (size_t iy = 0 ; iy < ny ; ++iy) {
      for (size_t i_spin = 0 ; i_spin < 2 ; ++i_spin) {
        for (size_t i_hop = 0 ; i_hop < n_hop ; ++i_hop) {
          auto v = vs[i_hop];
          auto dx = dxs[i_hop];
          auto dy = dys[i_hop];

          size_t jx = (nx + ix + dx) % nx;
          size_t jy = (ny + iy + dy) % ny;

          hop.add( v * cre(ix, iy, i_spin) * ann(jx, jy, i_spin));
        }
      }
    }
  }

  SectorGenerator<RepSize, SiteSize, Charge, Spin> sector_gen(system);
  auto max_qn = system.max_quantum_number();
  auto min_qn = system.min_quantum_number();

  std::ofstream outfile("basis.txt");
  for (auto charge = std::get<0>(min_qn); charge <= std::get<0>(max_qn); ++charge) {
    for (auto spin = std::get<1>(min_qn); spin <= std::get<1>(max_qn); ++spin) {
      auto sector = sector_gen.generate(charge, spin);
      size_t n_basis = sector.basis.size();
      if (n_basis ==  0) { continue; }

      cout << charge << "," << spin << endl;
      outfile << charge << "," << spin << endl;
      outfile << "=======================" << endl;
      for (auto v: sector.basis) { outfile << v << endl; }
      outfile << endl;

      Eigen::SparseMatrix<double> hamiltonian_matrix(n_basis, n_basis);
      {
        std::vector<Eigen::Triplet<double> > coefficients;

        for (size_t i_basis = 0; i_basis < n_basis; ++i_basis) {
          auto const &bvec_fvec = sector.basis[i_basis];
          auto row = hop.apply(std::get<0>(bvec_fvec), std::get<1>(bvec_fvec));
          for (auto const &r : row) {
            auto match_iter = sector.basismap.find(std::get<0>(r));
            if (match_iter == sector.basismap.end()) {
              cout << "ERROR!" << endl;
            } else {
              coefficients.emplace_back(match_iter->second, i_basis, std::get<2>(r));
            }
          } // for r in row (all resulting rows after applying H to col.
        } // for i_basis
        hamiltonian_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
      }

      cout << endl;
      cout << hamiltonian_matrix;
      cout << endl;
      cout << endl;
    }
  }
  outfile.close();
  return 0;

#if 0
  vector<tuple<bitset<RepSize>, bitset<SiteSize>>> basis_list;
  unordered_map<bitset<RepSize>, size_t> basis_map;
  for (auto iter = system.cbegin<RepSize, SiteSize>(Charge(1), Spin(1)) ; iter.valid() ; ++iter) {
    basis_list.push_back(iter.get());
    //cout << "BASIS GENERATED : " << iter.get() << endl;
  }

  for (size_t i = 0 ; i < basis_list.size() ; ++i) {
    basis_map[std::get<0>(basis_list[i])] = i;
    //cout << i << "\t" << std::get<0>(basis_list[i]) << endl;
  }


  std::vector<size_t> rows, cols;
  std::vector<double> vals;

  for (size_t i_basis = 0 ; i_basis < basis_list.size() ; ++i_basis) {
    auto const & bvec_fvec = basis_list[i_basis];
    auto row = hop.apply(std::get<0>(bvec_fvec), std::get<1>(bvec_fvec));
    for (auto const & r : row) {
      auto match_iter = basis_map.find(std::get<0>(r));
      //cout << std::get<0>(r) << "\t" << std::get<0>(bvec_fvec) << "\t" << std::get<2>(r) << endl;
      if (match_iter == basis_map.end()) {
        cout << "ERROR!" << endl;
      } else {

        rows.push_back(match_iter->second);
        cols.push_back(i_basis);
        vals.push_back(std::get<2>(r));
        //cout << match_iter->second << "\t" << i_basis << "\t" << std::get<2>(r) << endl;

        /*
        cout << "{" << match_iter->second;
        cout << "," << i_basis;
        cout << "," << std::get<2>(r);
        cout << "}," << endl;
         */
      }
    }
  }

  Eigen::MatrixXd mat(basis_list.size(), basis_list.size());

  mat.setZero();
  for (size_t i = 0 ; i < vals.size() ; ++i) {
  //  cout << rows[i] << "\t" << cols[i] << "\t" << vals[i] << endl;
    mat(rows[i], cols[i]) += vals[i];
  }
  Eigen::IOFormat heavy_fmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

  cout << mat.format(heavy_fmt) << endl;

#endif

  return 0;
}