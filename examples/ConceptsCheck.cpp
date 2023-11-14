#include <PlanetaryModel/All>
#include <concepts>
#include <iostream>
#include "PREM.h"

class DummyClass {};

template <typename T>
requires PlanetaryModel::SphericalGeometryModel<T>
void Dummy(T t) { t.UpperRadius(2); }

int main() {
    auto myprem = PREM();
    for (int i = 0; i < 13; ++i){
std::cout << myprem.Density(i)(0) << " " << myprem.VP(i)(0) << std::endl;
    }
    

}
