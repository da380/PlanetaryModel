#include <PlanetaryModel/All>
#include <concepts>
#include <iostream>

class DummyClass {};

template <typename T>
requires PlanetaryModel::SphericalGeometryModel<T>
void Dummy(T t) { t.UpperRadius(2); }

int main() {}
