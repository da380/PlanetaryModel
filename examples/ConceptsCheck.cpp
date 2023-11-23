#include <PlanetaryModel/All>
#include <concepts>
#include <iostream>

#include "EARTHMODELS.h"

class DummyClass {};

template <typename T>
    requires PlanetaryModel::SphericalGeometryModel<T>
void
Dummy(T t) {
    t.UpperRadius(2);
}

int
main() {
    using namespace EarthModels;
    auto myprem = EarthModels::PREM();
    auto mypertprem = EarthModels::PERTPREM();
    // //     for (int i = 0; i < 13; ++i){
    // // std::cout << myprem.Density(i)(0) << " " << myprem.VPH(i)(0) << " " <<
    // myprem.A(i)(0) << std::endl;
    // //     }
    //     std::cout << myprem.DensityNorm() << " "  << std::endl;
    static_assert(
        PlanetaryModel::HasNormalisationInformation<EarthConstants<double>>);
    static_assert(PlanetaryModel::SphericalGeometryModel<PREM<double, int>>);
    static_assert(PlanetaryModel::SphericalDensityModel<PREM<double, int>>);
    static_assert(PlanetaryModel::SphericalElasticModel<PREM<double, int>>);
    static_assert(
        PlanetaryModel::AsphericalDensityModel<PERTPREM<double, int>>);

    std::cout << myprem.VSH(2)(3480.0 / 6371.0) << std::endl;
    // polynomials:
}
