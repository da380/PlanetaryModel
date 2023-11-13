#ifndef PREM_MODEL_GUARD_H
#define PREM_MODEL_GUARD_H

#include <PlanetaryModel/All>
#include <cmath>
#include <concepts>
#include <functional>
#include <iostream>
#include <ranges>
#include <vector>

class EarthConstants {
 public:
  double LengthNorm() { return length_norm; };
  double MassNorm() { return mass_norm; };
  double TimeNorm() { return time_norm; }

  double DensityNorm() { return density_norm; };
  double InertiaNorm() { return inertia_norm; };
  double VelocityNorm() { return velocity_norm; };
  double AccelerationNorm() { return acceleration_norm; };
  double ForceNorm() { return force_norm; };
  double StressNorm() { return stress_norm; };
  double GravitationalConstant() { return gravitational_constant; };

 private:
  const double length_norm = 6.371 * std::pow(10.0, 6.0);
  const double mass_norm = 5.972 * std::pow(10.0, 24.0);
  const double time_norm = 3600.0;

  const double density_norm = mass_norm / std::pow(length_norm, 3.0);
  const double inertia_norm = mass_norm * std::pow(length_norm, 2.0);
  const double velocity_norm = length_norm / time_norm;
  const double acceleration_norm = length_norm / std::pow(time_norm, 2.0);
  const double force_norm = mass_norm * length_norm / std::pow(time_norm, 2.0);
  const double stress_norm =
      mass_norm / (std::pow(time_norm, 2.0) * length_norm);
  const double gravitational_constant =
      std::pow(length_norm, 3.0) / (mass_norm * std::pow(time_norm, 2.0));
};

class PREM {
 public:
  using size_type = int;
  using value_type = double;

  size_type NumberOfLayers() { return 13; };
  value_type LowerRadius(size_type i) { return vec_radii[i]; }
  value_type UpperRadius(size_type i) { return vec_radii[i + 1]; }
  value_type OuterRadius() { return vec_radii[13]; }

  // auto Density(size_type i) {return Density(i, value_type r);};
  value_type Density(size_type i, value_type r) {
    if (i == 0) {
      return 13.0885 - 8.8381 * std::pow(r / vec_radii[13], 2.0);
    }
  };

 private:
  std::vector<value_type> vec_radii{0.0,       1221500.0, 3480000.0, 3630000.0,
                                    5600000.0, 5701000.0, 5771000.0, 5971000.0,
                                    6151000.0, 6291000.0, 6346600.0, 6356000.0,
                                    6368000.0, 6371000.0};
};

#endif