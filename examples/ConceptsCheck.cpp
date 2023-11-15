#include <PlanetaryModel/All>
#include <concepts>
#include <iostream>

#include "PREM.h"

class DummyClass {};

template <typename T>
requires PlanetaryModel::SphericalGeometryModel<T>
void Dummy(T t) {
  t.UpperRadius(2);
}

int main() {
  auto myprem = PREM();
  // //     for (int i = 0; i < 13; ++i){
  // // std::cout << myprem.Density(i)(0) << " " << myprem.VPH(i)(0) << " " <<
  // myprem.A(i)(0) << std::endl;
  // //     }
  //     std::cout << myprem.DensityNorm() << " "  << std::endl;
  static_assert(
      PlanetaryModel::HasNormalisationInformation<EarthConstants<double> >);
  static_assert(PlanetaryModel::SphericalGeometryModel<PREM<double, int> >);
  static_assert(PlanetaryModel::SphericalDensityModel<PREM<double, int> >);
  static_assert(PlanetaryModel::SphericalElasticModel<PREM<double, int> >);

  // polynomials:
  Interpolation::Polynomial1D<double> vecpoly{1.0, 2.0};
  Interpolation::Polynomial1D<double> vecpoly2 = vecpoly * 2.0;
  Interpolation::Polynomial1D<double> vecpoly3 = vecpoly + 2.0;
  Interpolation::Polynomial1D<double> vecpoly4 = vecpoly + vecpoly2;
  std::cout << vecpoly(0) << " " << vecpoly2(0) << std::endl;
  std::cout << vecpoly(0) << " " << vecpoly3(0) << std::endl;
  std::cout << vecpoly(0) << "," << vecpoly(1) << " " << vecpoly3(0) << "," << vecpoly3(1) << std::endl;
  //   vecpoly2 = vecpoly * 2.0;
}
