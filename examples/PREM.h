#ifndef PREM_MODEL_GUARD_H
#define PREM_MODEL_GUARD_H

#include "Interpolation/All"
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
    const double force_norm =
        mass_norm * length_norm / std::pow(time_norm, 2.0);
    const double stress_norm =
        mass_norm / (std::pow(time_norm, 2.0) * length_norm);
    const double gravitational_constant =
        std::pow(length_norm, 3.0) / (mass_norm * std::pow(time_norm, 2.0));
};

class PREM {
  public:
    using size_type = int;
    using value_type = double;

    // Constructor
    PREM(){};

    // functions
    size_type NumberOfLayers() { return 13; };
    value_type LowerRadius(size_type i) { return vec_radii[i]; }
    value_type UpperRadius(size_type i) { return vec_radii[i + 1]; }
    value_type OuterRadius() { return vec_radii[13]; }

    Interpolation::Polynomial1D<double> Density(size_type i) {
        return vec_density[i];
    };
    Interpolation::Polynomial1D<double> VP(size_type i) {
        return vec_p_velocity[i];
    };
    Interpolation::Polynomial1D<double> VPV(size_type i) {
        return vec_pv_velocity[i];
    };
    Interpolation::Polynomial1D<double> VPH(size_type i) {
        return vec_ph_velocity[i];
    };
    Interpolation::Polynomial1D<double> VS(size_type i) {
        return vec_p_velocity[i];
    };
    Interpolation::Polynomial1D<double> VSV(size_type i) {
        return vec_pv_velocity[i];
    };
    Interpolation::Polynomial1D<double> VSH(size_type i) {
        return vec_ph_velocity[i];
    };

    // data
  private:
    std::vector<value_type> vec_radii{
        0.0,       1221500.0, 3480000.0, 3630000.0, 5600000.0,
        5701000.0, 5771000.0, 5971000.0, 6151000.0, 6291000.0,
        6346600.0, 6356000.0, 6368000.0, 6371000.0};

    std::vector<Interpolation::Polynomial1D<double>> vec_density{
        {13.0855, 0, -8.8381},
        {12.5815, -1.2638, -3.6426, -5.5281},
        {7.9565, -6.4761, 5.5283, -3.0807},
        {7.9565, -6.4761, 5.5283, -3.0807},
        {7.9565, -6.4761, 5.5283, -3.0807},
        {5.3197, -1.4836},
        {11.2494, -8.0298},
        {7.1089, -3.8045},
        {2.6910, 0.6924},
        {2.6910, 0.6924},
        {2.900},
        {2.600},
        {1.020}};

    std::vector<Interpolation::Polynomial1D<double>> vec_p_velocity{
        {11.2622, 0, -6.3640},
        {11.0487, -4.0362, 4.8023, -13.5732},
        {15.3891, -5.3181, 5.5242, -2.5514},
        {24.9520, -40.4673, 51.4832, -26.6419},
        {29.2766, -23.6027, 5.5242, -2.5514},
        {19.0957, -9.8672},
        {39.7027, -32.6166},
        {20.3926, -12.2569},
        {4.1875, 3.9382},
        {4.1875, 3.9382},
        {6.800},
        {5.800},
        {1.450}};
    std::vector<Interpolation::Polynomial1D<double>> vec_pv_velocity{
        {11.2622, 0, -6.3640},
        {11.0487, -4.0362, 4.8023, -13.5732},
        {15.3891, -5.3181, 5.5242, -2.5514},
        {24.9520, -40.4673, 51.4832, -26.6419},
        {29.2766, -23.6027, 5.5242, -2.5514},
        {19.0957, -9.8672},
        {39.7027, -32.6166},
        {20.3926, -12.2569},
        {0.8317, 7.2180},
        {0.8317, 7.2180},
        {6.800},
        {5.800},
        {1.450}};
    std::vector<Interpolation::Polynomial1D<double>> vec_ph_velocity{
        {11.2622, 0, -6.3640},
        {11.0487, -4.0362, 4.8023, -13.5732},
        {15.3891, -5.3181, 5.5242, -2.5514},
        {24.9520, -40.4673, 51.4832, -26.6419},
        {29.2766, -23.6027, 5.5242, -2.5514},
        {19.0957, -9.8672},
        {39.7027, -32.6166},
        {20.3926, -12.2569},
        {3.5908, 4.6172},
        {3.5908, 4.6172},
        {6.800},
        {5.800},
        {1.450}};

    std::vector<Interpolation::Polynomial1D<double>> vec_s_velocity{
        {3.6678, 0, -4.4475},
        {0.0},
        {6.9254, 1.4672, -2.0834, 0.9783},
        {11.1671, -13.7818, 17.4575, -9.2777},
        {22.3459, -17.2473, -2.0834, 0.9783},
        {9.9839, -4.9324},
        {22.3512, -18.5856},
        {8.9496, -4.4597},
        {2.1519, 2.3481},
        {2.1519, 2.3481},
        {3.900},
        {3.200},
        {0}};
    std::vector<Interpolation::Polynomial1D<double>> vec_sv_velocity{
        {3.6678, 0, -4.4475},
        {0.0},
        {6.9254, 1.4672, -2.0834, 0.9783},
        {11.1671, -13.7818, 17.4575, -9.2777},
        {22.3459, -17.2473, -2.0834, 0.9783},
        {9.9839, -4.9324},
        {22.3512, -18.5856},
        {8.9496, -4.4597},
        {5.8582, -1.4678},
        {5.8582, -1.4678},
        {3.900},
        {3.200},
        {0}};
    std::vector<Interpolation::Polynomial1D<double>> vec_sh_velocity{
        {3.6678, 0, -4.4475},
        {0.0},
        {6.9254, 1.4672, -2.0834, 0.9783},
        {11.1671, -13.7818, 17.4575, -9.2777},
        {22.3459, -17.2473, -2.0834, 0.9783},
        {9.9839, -4.9324},
        {22.3512, -18.5856},
        {8.9496, -4.4597},
        {-1.0839, 5.7176},
        {-1.0839, 5.7176},
        {3.900},
        {3.200},
        {0}};
    std::vector<Interpolation::Polynomial1D<double>> vec_Qmu{
        {84.6}, {std::pow(10.0, 10.0)},
        {312},  {312},
        {312},  {143},
        {143},  {143},
        {80},   {600},
        {600},  {std::pow(10.0, 10.0)}};

    std::vector<Interpolation::Polynomial1D<double>> vec_QKappa{
        {1327.7}, {57823}, {57823}, {57823}, {57823}, {57823},
        {57823},  {57823}, {57823}, {57823}, {57823}, {57823}};
};

#endif