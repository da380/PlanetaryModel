#ifndef EARTHMODELS_GUARD_H
#define EARTHMODELS_GUARD_H

// #include <PlanetaryModel/All>
// #include <cmath>
// #include <concepts>
#include <fstream>
// #include <functional>
// #include <iostream>
// #include <ranges>
// #include <vector>

// #include "Interpolation/All"
#include <Interpolation/CubicSpline>
#include <Interpolation/Polynomial>

namespace EarthModels {
template <typename FLOAT = double> class EarthConstants {
  public:
    using value_type = FLOAT;
    FLOAT LengthNorm() const { return length_norm; };
    FLOAT MassNorm() const { return mass_norm; };
    FLOAT TimeNorm() const { return time_norm; }

    FLOAT DensityNorm() const { return density_norm; };
    FLOAT InertiaNorm() const { return inertia_norm; };
    FLOAT VelocityNorm() const { return velocity_norm; };
    FLOAT AccelerationNorm() const { return acceleration_norm; };
    FLOAT ForceNorm() const { return force_norm; };
    FLOAT StressNorm() const { return stress_norm; };
    FLOAT GravitationalConstant() const { return gravitational_constant; };

  private:
    const FLOAT length_norm = 6.371 * std::pow(10.0, 6.0);
    const FLOAT mass_norm = 5.972 * std::pow(10.0, 24.0);
    const FLOAT time_norm = 3600.0;

    const FLOAT density_norm = mass_norm / std::pow(length_norm, 3.0);
    const FLOAT inertia_norm = mass_norm * std::pow(length_norm, 2.0);
    const FLOAT velocity_norm = length_norm / time_norm;
    const FLOAT acceleration_norm = length_norm / std::pow(time_norm, 2.0);
    const FLOAT force_norm = mass_norm * length_norm / std::pow(time_norm, 2.0);
    const FLOAT stress_norm =
        mass_norm / (std::pow(time_norm, 2.0) * length_norm);
    const FLOAT gravitational_constant =
        std::pow(length_norm, 3.0) / (mass_norm * std::pow(time_norm, 2.0));
};

template <typename FLOAT = double, typename INTEGRAL = int>
class PREM : public EarthConstants<FLOAT> {
  public:
    using size_type = INTEGRAL;

    // Constructor
    PREM() {
        for (int idx = 0; idx < 13; ++idx) {
            vec_A.push_back(vec_density[idx] * vec_ph_velocity[idx] *
                            vec_ph_velocity[idx]);
        };
    };

    // Geometry of PREM
    INTEGRAL NumberOfLayers() const { return 13; };
    FLOAT LowerRadius(INTEGRAL i) const {
        return vec_radii[i] / this->LengthNorm();
    }
    FLOAT UpperRadius(INTEGRAL i) const {
        return vec_radii[i + 1] / this->LengthNorm();
    }
    FLOAT OuterRadius() const { return vec_radii[13] / this->LengthNorm(); }

    // Density
    Interpolation::Polynomial1D<FLOAT> Density(INTEGRAL i) const {
        if (i < 0) {
            Interpolation::Polynomial1D<FLOAT> rettemp{0.0};
            return rettemp;
        } else if (i > this->NumberOfLayers() - 1) {
            Interpolation::Polynomial1D<FLOAT> rettemp{0.0};
            return rettemp;
        } else {
            return 1000.0 * vec_density[i] / this->DensityNorm();
        }
    };

    // Isotropy/fluid/solid etc
    bool IsIsotropic() const { return false; };

    // Solid or fluid
    bool IsSolid(INTEGRAL i) const {
        if (i == 1 || i == 12) {
            return false;
        } else {
            return true;
        }
    }
    bool IsFluid(INTEGRAL i) const { return !IsSolid(i); }

    // Return TI elastic modulii

    // Velocities

    Interpolation::Polynomial1D<FLOAT> VP(INTEGRAL i) const {
        return vec_p_velocity[i] * 1000.0 / this->VelocityNorm();
    };
    Interpolation::Polynomial1D<FLOAT> VPV(INTEGRAL i) const {
        return vec_pv_velocity[i] * 1000.0 / this->VelocityNorm();
    };
    Interpolation::Polynomial1D<FLOAT> VPH(INTEGRAL i) const {
        return vec_ph_velocity[i] * 1000.0 / this->VelocityNorm();
    };
    Interpolation::Polynomial1D<FLOAT> VS(INTEGRAL i) const {
        return vec_s_velocity[i] * 1000.0 / this->VelocityNorm();
    };
    Interpolation::Polynomial1D<FLOAT> VSV(INTEGRAL i) const {
        return vec_sv_velocity[i] * 1000.0 / this->VelocityNorm();
    };
    Interpolation::Polynomial1D<FLOAT> VSH(INTEGRAL i) const {
        return vec_sh_velocity[i] * 1000.0 / this->VelocityNorm();
    };

    // Returning eta, A, C, N, L, kappa, mu
    auto Eta(INTEGRAL i) const { return vec_eta[i]; }
    auto A(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VPH(i)(x) * VPH(i)(x);
        };
        // auto aret = Density(i) * VPH(i) * VPH(i);
        return aret;
    };
    auto C(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VPV(i)(x) * VPV(i)(x);
        };
        return aret;
    };
    auto N(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VSH(i)(x) * VSH(i)(x);
        };
        return aret;
    };
    auto L(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VSV(i)(x) * VSV(i)(x);
        };
        return aret;
    };
    auto F(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Eta(i)(x) * (A(i)(x) - 2 * L(i)(x));
        };
        return aret;
    };
    auto Kappa(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return (C(i)(x) + 4.0 * (A(i)(x) - N(i)(x) + F(i)(x))) / 9.0;
        };
        return aret;
    };
    auto Mu(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return (C(i)(x) + A(i)(x) + 6.0 * L(i)(x) + 5.0 * N(i)(x) -
                    2.0 * F(i)(x)) /
                   15.0;
        };
        return aret;
    };

    // data
  private:
    std::vector<FLOAT> vec_radii{0.0,       1221500.0, 3480000.0, 3630000.0,
                                 5600000.0, 5701000.0, 5771000.0, 5971000.0,
                                 6151000.0, 6291000.0, 6346600.0, 6356000.0,
                                 6368000.0, 6371000.0};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_density{
        {13.0885, 0, -8.8381},
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

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_p_velocity{
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
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pv_velocity{
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
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_ph_velocity{
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

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_s_velocity{
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
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_sv_velocity{
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
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_sh_velocity{
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
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_Qmu{
        {84.6}, {std::pow(10.0, 10.0)},
        {312},  {312},
        {312},  {143},
        {143},  {143},
        {80},   {600},
        {600},  {std::pow(10.0, 10.0)}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_QKappa{
        {1327.7}, {57823}, {57823}, {57823}, {57823}, {57823},
        {57823},  {57823}, {57823}, {57823}, {57823}, {57823}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_eta{{1},
                                                            {1},
                                                            {1},
                                                            {1},
                                                            {1},
                                                            {1},
                                                            {1},
                                                            {1},
                                                            {1},
                                                            {1},
                                                            {3.3687, -2.4778},
                                                            {3.3687, -2.4778},
                                                            {1},
                                                            {1},
                                                            {1}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_A;
};

template <typename FLOAT = double, typename INTEGRAL = int>
class PERTPREM : public PREM<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    PERTPREM() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const { return 0.0; };
    FLOAT MaxRadius() const { return this->OuterRadius(); };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0},
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0},
        {0, 0, 0}, {0, 0, 0}, {0.0}};
};

template <typename FLOAT> class HomogeneousConstants {
  public:
    using value_type = FLOAT;
    FLOAT LengthNorm() const { return length_norm; };
    FLOAT MassNorm() const { return mass_norm; };
    FLOAT TimeNorm() const { return time_norm; }

    FLOAT DensityNorm() const { return density_norm; };
    FLOAT InertiaNorm() const { return inertia_norm; };
    FLOAT VelocityNorm() const { return velocity_norm; };
    FLOAT AccelerationNorm() const { return acceleration_norm; };
    FLOAT ForceNorm() const { return force_norm; };
    FLOAT StressNorm() const { return stress_norm; };
    FLOAT GravitationalConstant() const { return gravitational_constant; };

  private:
    const FLOAT length_norm = 6.371 * std::pow(10.0, 6.0);
    const FLOAT mass_norm = 5.972 * std::pow(10.0, 24.0);
    const FLOAT time_norm = 3600.0;

    const FLOAT density_norm = mass_norm / std::pow(length_norm, 3.0);
    const FLOAT inertia_norm = mass_norm * std::pow(length_norm, 2.0);
    const FLOAT velocity_norm = length_norm / time_norm;
    const FLOAT acceleration_norm = length_norm / std::pow(time_norm, 2.0);
    const FLOAT force_norm = mass_norm * length_norm / std::pow(time_norm, 2.0);
    const FLOAT stress_norm =
        mass_norm / (std::pow(time_norm, 2.0) * length_norm);
    const FLOAT gravitational_constant =
        std::pow(length_norm, 3.0) / (mass_norm * std::pow(time_norm, 2.0));
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOSPHERE : public HomogeneousConstants<FLOAT> {
  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOSPHERE() {
        for (int idx = 0; idx < 3; ++idx) {
            vec_A.push_back(vec_density[idx] * vec_ph_velocity[idx] *
                            vec_ph_velocity[idx]);
        };
    };

    // Geometry of PREM
    INTEGRAL NumberOfLayers() const { return 3; };
    FLOAT LowerRadius(INTEGRAL i) const {
        return vec_radii[i] / this->LengthNorm();
    }
    FLOAT UpperRadius(INTEGRAL i) const {
        return vec_radii[i + 1] / this->LengthNorm();
    }
    FLOAT OuterRadius() const { return vec_radii[3] / this->LengthNorm(); }

    // Density
    Interpolation::Polynomial1D<FLOAT> Density(INTEGRAL i) const {
        if (i < 0) {
            Interpolation::Polynomial1D<FLOAT> rettemp{0.0};
            return rettemp;
        } else if (i > this->NumberOfLayers() - 1) {
            Interpolation::Polynomial1D<FLOAT> rettemp{0.0};
            return rettemp;
        } else {
            return 1000.0 * vec_density[i] / this->DensityNorm();
        }
    };

    // Isotropy/fluid/solid etc
    bool IsIsotropic() const { return false; };

    // Solid or fluid
    bool IsSolid(INTEGRAL i) const { return true; }
    bool IsFluid(INTEGRAL i) const { return !IsSolid(i); }

    // Return TI elastic modulii

    // Velocities

    Interpolation::Polynomial1D<FLOAT> VP(INTEGRAL i) const {
        return vec_p_velocity[i];
    };
    Interpolation::Polynomial1D<FLOAT> VPV(INTEGRAL i) const {
        return vec_pv_velocity[i];
    };
    Interpolation::Polynomial1D<FLOAT> VPH(INTEGRAL i) const {
        return vec_ph_velocity[i];
    };
    Interpolation::Polynomial1D<FLOAT> VS(INTEGRAL i) const {
        return vec_s_velocity[i];
    };
    Interpolation::Polynomial1D<FLOAT> VSV(INTEGRAL i) const {
        return vec_sv_velocity[i];
    };
    Interpolation::Polynomial1D<FLOAT> VSH(INTEGRAL i) const {
        return vec_sh_velocity[i];
    };

    // Returning eta, A, C, N, L, kappa, mu
    auto Eta(INTEGRAL i) const { return vec_eta[i]; }
    auto A(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VPH(i)(x) * VPH(i)(x);
        };
        // auto aret = Density(i) * VPH(i) * VPH(i);
        return aret;
    };
    auto C(INTEGRAL i) const {
        // auto aret = [i, this](FLOAT x) {
        //     return Density(i)(x) * VPV(i)(x) * VPV(i)(x);
        // };
        return vec_A[i];
    };
    auto N(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VSH(i)(x) * VSH(i)(x);
        };
        return aret;
    };
    auto L(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VSV(i)(x) * VSV(i)(x);
        };
        return aret;
    };
    auto F(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Eta(i)(x) * (A(i)(x) - 2 * L(i)(x));
        };
        return aret;
    };
    auto Kappa(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return (C(i)(x) + 4.0 * (A(i)(x) - N(i)(x) + F(i)(x))) / 9.0;
        };
        return aret;
    };
    auto Mu(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return (C(i)(x) + A(i)(x) + 6.0 * L(i)(x) + 5.0 * N(i)(x) -
                    2.0 * F(i)(x)) /
                   15.0;
        };
        return aret;
    };

    // data
  private:
    std::vector<FLOAT> vec_radii{0.0, 1221500.0, 2700000.0, 6371000.0};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_density{
        {5.51, 0, 0.0}, {5.51, 0, 0.0}, {5.51, 0, 0.0}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_p_velocity{
        {1.0, 0, 0}, {1.0, 0, 0}, {1.0, 0, 0}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pv_velocity{
        {1.0, 0, 0}, {1.0, 0, 0}, {1.0, 0, 0}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_ph_velocity{
        {1.0, 0, 0}, {1.0, 0, 0}, {1.0, 0, 0}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_s_velocity{
        {1.0, 0, 0}, {1.0, 0, 0}, {1.0, 0, 0}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_sv_velocity{
        {1.0, 0, 0}, {1.0, 0, 0}, {1.0, 0, 0}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_sh_velocity{
        {1.0, 0, 0}, {1.0, 0, 0}, {1.0, 0, 0}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_Qmu{
        {100.0}, {100.0}, {100.0}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_QKappa{
        {57823}, {57823}, {57823}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_eta{{1}, {1}, {1}};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_A;
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND0 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND0() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const { return 0.0; };
    FLOAT MaxRadius() const { return this->OuterRadius(); };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND1 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND1() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const { return 0.01 * r; };
    FLOAT MaxRadius() const { return 1.01 * this->OuterRadius(); };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND2 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND2() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const { return 0.02 * r; };
    FLOAT MaxRadius() const { return 1.02 * this->OuterRadius(); };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND3 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND3() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const {
        return 0.2 * this->OuterRadius() * (r / this->OuterRadius()) *
               (1.0 - r / this->OuterRadius());
    };
    FLOAT MaxRadius() const { return 1.0 * this->OuterRadius(); };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND4 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND4() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const {
        return 0.8 * this->OuterRadius() * (r / this->OuterRadius()) *
               (1.0 - r / this->OuterRadius());
    };
    FLOAT MaxRadius() const { return 1.0 * this->OuterRadius(); };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND5 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND5() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const {
        return r * r - r;
    };
    FLOAT MaxRadius() const {
        return this->OuterRadius() * this->OuterRadius();
    };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND6 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND6() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const {
        return 0.01 * r * (1 - r / this->OuterRadius()) * std::sin(theta);
    };
    FLOAT MaxRadius() const { return this->OuterRadius() * 1.0; };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND7 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND7() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const {
        return 0.5 * r * (1 - r / this->OuterRadius()) * std::sin(theta);
    };
    FLOAT MaxRadius() const { return this->OuterRadius() * 1.0; };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND8 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND8() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const {
        return 0.5 * r * (1 - r / this->OuterRadius()) * std::cos(theta);
    };
    FLOAT MaxRadius() const { return this->OuterRadius() * 1.0; };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND9 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND9() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const {
        return 0.5 * r * (1 - r / this->OuterRadius()) * std::sin(theta) *
               std::sin(phi);
    };
    FLOAT MaxRadius() const { return this->OuterRadius() * 1.0; };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int>
class HOMOBOUND10 : public HOMOSPHERE<FLOAT, int> {

  public:
    using size_type = INTEGRAL;

    // Constructor
    HOMOBOUND10() {};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const {
        return 0.02 * r * (1 - r / this->OuterRadius()) * std::cos(theta);
    };
    FLOAT MaxRadius() const { return this->OuterRadius() * 1.0; };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

template <typename FLOAT = double, typename INTEGRAL = int> class ModelInput {
  public:
    using size_type = INTEGRAL;
    using value_type = FLOAT;
    // using InterpA = Interpolation::Akima<std::vector<double>::iterator,
    //                                      std::vector<double>::iterator>;
    // template <typename NFLOAT = FLOAT>
    // using InterpA = Interpolation::CubicSpline<std::vector<double>::iterator,
    //                                            std::vector<double>::iterator>;
    using myvector = std::vector<FLOAT>;
    using myiter = myvector::iterator;
    using InterpA = Interpolation::CubicSpline<myiter, myiter>;
    // using InterpA = Interpolation::CubicSpline<std::vector<FLOAT>::iterator,
    //                                            std::vector<FLOAT>::iterator>;
    // constructors
    ModelInput() {};
    ModelInput(const std::string &);
    template <template <typename> class ParameterModel>
    ModelInput(const std::string &, const ParameterModel<FLOAT> &);

    // norms
    FLOAT LengthNorm() const { return length_norm; };
    FLOAT MassNorm() const { return mass_norm; };
    FLOAT TimeNorm() const { return time_norm; }
    FLOAT DensityNorm() const { return density_norm; };
    FLOAT InertiaNorm() const { return inertia_norm; };
    FLOAT VelocityNorm() const { return velocity_norm; };
    FLOAT AccelerationNorm() const { return acceleration_norm; };
    FLOAT ForceNorm() const { return force_norm; };
    FLOAT StressNorm() const { return stress_norm; };
    FLOAT GravitationalConstant() const { return gravitational_constant; };

    // Geometry of model
    int NumberOfLayers() const { return _numlayers; };
    int LayerLowerIndex(int i) const { return _vec_indices[i][0]; };
    int LayerUpperIndex(int i) const { return _vec_indices[i][1]; };
    int LayerIndexDifference(int i) const {
        return _vec_indices[i][1] - _vec_indices[i][0];
    };
    auto LayerRadii() const { return layered_radii; };
    auto LayerRadii(int i) const { return layered_radii[i]; };
    auto LayerRadiiNumber(int i) const { return layered_radii[i].size(); };

    FLOAT LowerRadius(INTEGRAL i) const {
        if (i < 0) {
            throw std::invalid_argument("Negative layer index");
        } else if (i > _numlayers - 1) {
            assert("Outside the number of layers in the model");
            throw std::invalid_argument(
                "Layer index greater than number of layers");
        };
        return vec_layers[i];
    }
    FLOAT UpperRadius(INTEGRAL i) const {
        if (i < 0) {
            throw std::invalid_argument("Negative layer index");
        } else if (i > _numlayers - 1) {
            assert("Outside the number of layers in the model");
            throw std::invalid_argument(
                "Layer index greater than number of layers");
        };
        return vec_layers[i + 1];
    }
    FLOAT OuterRadius() const { return vec_layers[_numlayers]; }

    // Isotropy/fluid/solid etc
    bool IsIsotropic() const { return _isisotropic; };

    // Solid or fluid
    bool IsSolid(INTEGRAL i) const { return _issolid[i]; }
    bool IsFluid(INTEGRAL i) const { return !IsSolid(i); }

    // Density
    InterpA Density(INTEGRAL i) const {
        if (i < 0) {
            throw std::invalid_argument("Negative layer index");
        } else if (i > _numlayers - 1) {
            assert("Outside model");
            throw std::invalid_argument(
                "Layer index greater than number of layers");
        };
        return func_rho[i];
        // return func_rhoc[i];
    };

    // Velocities
    // InterpA VP(INTEGRAL i) const { return func_vpv[i]; };
    InterpA VPV(INTEGRAL i) const {
        if (i < 0) {
            throw std::invalid_argument("Negative layer index");
        } else if (i > _numlayers - 1) {
            assert("Outside model");
            throw std::invalid_argument(
                "Layer index greater than number of layers");
        };

        return func_vpv[i];
    };
    InterpA VPH(INTEGRAL i) const {
        if (i < 0) {
            throw std::invalid_argument("Negative layer index");
        } else if (i > _numlayers - 1) {
            assert("Outside model");
            throw std::invalid_argument(
                "Layer index greater than number of layers");
        };
        return func_vph[i];
    };
    // InterpA VS(INTEGRAL i) const { return vec_s_velocity[i]; };
    InterpA VSV(INTEGRAL i) const {
        if (i < 0) {
            throw std::invalid_argument("Negative layer index");
        } else if (i > _numlayers - 1) {
            assert("Outside model");
            throw std::invalid_argument(
                "Layer index greater than number of layers");
        };
        return func_vsv[i];
    };
    InterpA VSH(INTEGRAL i) const {
        if (i < 0) {
            throw std::invalid_argument("Negative layer index");
        } else if (i > _numlayers - 1) {
            assert("Outside model");
            throw std::invalid_argument(
                "Layer index greater than number of layers");
        };
        return func_vsh[i];
    };

    ///////////////////////////////////////////////////////////////
    ////////////////// !!!!!!!!!!!!!!!!!!!!!!!!!!! ////////////////
    ///////////////////////////////////////////////////////////////
    InterpA VS(INTEGRAL i) const {
        if (i < 0) {
            throw std::invalid_argument("Negative layer index");
        } else if (i > _numlayers - 1) {
            assert("Outside model");
            throw std::invalid_argument(
                "Layer index greater than number of layers");
        };
        return func_vsv[i];
    };
    InterpA VP(INTEGRAL i) const {
        if (i < 0) {
            throw std::invalid_argument("Negative layer index");
        } else if (i > _numlayers - 1) {
            assert("Outside model");
            throw std::invalid_argument(
                "Layer index greater than number of layers");
        };
        return func_vpv[i];
    };
    ///////////////////////////////////////////////////////////////
    ////////////////// !!!!!!!!!!!!!!!!!!!!!!!!!!! ////////////////
    ///////////////////////////////////////////////////////////////

    // Returning eta:
    auto Eta(INTEGRAL i) const {
        if (i < 0) {
            throw std::invalid_argument("Negative layer index");
        } else if (i > _numlayers - 1) {
            assert("Outside model");
            throw std::invalid_argument(
                "Layer index greater than number of layers");
        };
        return func_eta[i];
    }

    // returning A, C, N, L, kappa, mu
    auto A(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VPH(i)(x) * VPH(i)(x);
        };
        return aret;
    };
    auto C(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VPV(i)(x) * VPV(i)(x);
        };
        return aret;
    };
    auto N(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VSH(i)(x) * VSH(i)(x);
        };
        return aret;
    };
    auto L(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Density(i)(x) * VSV(i)(x) * VSV(i)(x);
        };
        return aret;
    };
    auto F(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return Eta(i)(x) * (A(i)(x) - 2 * L(i)(x));
        };
        return aret;
    };
    auto Kappa(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return (C(i)(x) + 4.0 * (A(i)(x) - N(i)(x) + F(i)(x))) / 9.0;
        };
        return aret;
    };
    auto Mu(INTEGRAL i) const {
        auto aret = [i, this](FLOAT x) {
            return (C(i)(x) + A(i)(x) + 6.0 * L(i)(x) + 5.0 * N(i)(x) -
                    2.0 * F(i)(x)) /
                   15.0;
        };
        return aret;
    };

  private:
    using vecdb = std::vector<double>;
    using vvecdb = std::vector<vecdb>;
    vecdb vec_radius, vec_rho, vec_vpv, vec_vsv, vec_qkappa, vec_qshear,
        vec_vph, vec_vsh, vec_eta, vec_layers;
    vvecdb layered_radii = vvecdb(1, vecdb());
    vvecdb layered_rho = vvecdb(1, vecdb());
    vvecdb layered_vpv = vvecdb(1, vecdb());
    vvecdb layered_vsv = vvecdb(1, vecdb());
    vvecdb layered_qkappa = vvecdb(1, vecdb());
    vvecdb layered_qshear = vvecdb(1, vecdb());
    vvecdb layered_vph = vvecdb(1, vecdb());
    vvecdb layered_vsh = vvecdb(1, vecdb());
    vvecdb layered_eta = vvecdb(1, vecdb());

    bool _isisotropic = true;
    std::vector<bool> _issolid = std::vector<bool>(1, true);
    std::vector<std::vector<int>> _vec_indices;
    // density
    // Interpolation::CubicSpline<std::vector<double>::iterator,
    //                            std::vector<double>::iterator>
    //     checkval;
    std::vector<InterpA> func_rho, func_vpv, func_vsv, func_qkappa, func_qshear,
        func_vph, func_vsh, func_eta;
    // std::vector<InterpC> func_rhoc;
    // vvecdb layered_radii = vvecdb(1, vecdb());
    // std::vector<std::vector<double>> layered_vpv, layered_vsv,
    //     layered_qkappa, layered_qshear, layered_vph, layered_vsh,
    //     layered_eta;

    // model information
    std::string modeltitle;
    int ifanis, ifdeck, numnodes, nic, noc, _numlayers;
    double tref;

    FLOAT length_norm, mass_norm, time_norm, density_norm, inertia_norm,
        velocity_norm, acceleration_norm, force_norm, stress_norm,
        gravitational_constant;

    ///////////////////////////////////////////////
    /////////////// ??????????????? ///////////////
    ///////////////////////////////////////////////

    // find layers of model
    std::vector<double> findlayers(const std::vector<double> &vec_sorted) {
        std::vector<double> vec_bounds;
        auto i1 = vec_sorted.begin();
        while (i1 != vec_sorted.end()) {
            vec_bounds.push_back(*i1);
            i1 = std::adjacent_find(++i1, vec_sorted.end());
        }
        vec_bounds.push_back(*(--i1));
        return vec_bounds;
    };

    // find indices
    std::vector<std::size_t>
    layerindices(const std::vector<double> &vec_sorted) {
        std::vector<std::size_t> vec_indices;
        auto i1 = vec_sorted.begin();
        while (i1 != vec_sorted.end()) {
            vec_indices.push_back(std::distance(vec_sorted.begin(), i1));
            i1 = std::adjacent_find(++i1, vec_sorted.end());
        }
        vec_indices.push_back(
            std::distance(vec_sorted.begin(), vec_sorted.end()) - 1);
        return vec_indices;
    }
};

// "default" constructor
template <typename FLOAT, typename INTEGRAL>
ModelInput<FLOAT, INTEGRAL>::ModelInput(const std::string &pathtofile)
    : ModelInput(pathtofile, EarthModels::EarthConstants<FLOAT>()){};

// full constructor
template <typename FLOAT, typename INTEGRAL>
template <template <typename> class ParameterModel>
ModelInput<FLOAT, INTEGRAL>::ModelInput(
    const std::string &pathtofile, const ParameterModel<FLOAT> &ModelConstants)
    : length_norm(ModelConstants.LengthNorm()),
      mass_norm(ModelConstants.MassNorm()),
      time_norm(ModelConstants.TimeNorm()),
      density_norm(ModelConstants.DensityNorm()),
      velocity_norm(ModelConstants.VelocityNorm()),
      acceleration_norm(ModelConstants.AccelerationNorm()),
      force_norm(ModelConstants.ForceNorm()),
      stress_norm(ModelConstants.StressNorm()),
      inertia_norm(ModelConstants.InertiaNorm()),
      gravitational_constant(ModelConstants.GravitationalConstant()) {

    // std::cout << this->density_norm << "\n";
    // opening file
    std::fstream modelfile;
    modelfile.open(pathtofile, std::ios::in);

    // getting information out of file
    if (modelfile.is_open()) {
        // get first line (title)
        getline(modelfile, modeltitle);

        // extract information from second line and move to next line
        modelfile >> ifanis >> tref >> ifdeck;
        modelfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // extract information from third line and move to next line
        modelfile >> numnodes >> nic >> noc;
        modelfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // loop through the deck
        int laynum = 0;
        int idxouter = 0;
        int idxinner = 0;

        while (idxouter < numnodes) {

            std::vector<double> tmp_radius;
            bool samelayer = true;

            // while (samelayer) {
            // double radius;
            double radius, rho, vpv, vsv, qkappa, qshear, vph, vsh, eta;
            modelfile >> radius >> rho >> vpv >> vsv >> qkappa >> qshear >>
                vph >> vsh >> eta;
            if (idxinner > 0 && (radius / this->length_norm ==
                                 layered_radii[laynum][idxinner - 1])) {
                // move to next layer
                layered_radii.push_back({radius / this->length_norm});
                layered_rho.push_back({rho / this->density_norm});
                layered_vpv.push_back({vpv / this->velocity_norm});
                layered_vsv.push_back({vsv / this->velocity_norm});
                layered_qkappa.push_back({qkappa});
                layered_qshear.push_back({qshear});
                layered_vph.push_back({vph / this->velocity_norm});
                layered_vsh.push_back({vsh / this->velocity_norm});
                layered_eta.push_back({eta});

                // isotropy
                if (_isisotropic && vpv != vsv) {
                    _isisotropic = false;
                }

                // fluid/solid:
                if (vsv == 0.0 && vsh == 0) {
                    _issolid.push_back(false);
                } else {
                    _issolid.push_back(true);
                }

                // set idxinner back to zero
                idxinner = 0;
                ++laynum;
            } else {
                // put next value in current layer in
                layered_radii[laynum].push_back(radius / this->length_norm);
                layered_rho[laynum].push_back(rho / this->density_norm);
                layered_vpv[laynum].push_back(vpv / this->velocity_norm);
                layered_vsv[laynum].push_back(vsv / this->velocity_norm);
                layered_qkappa[laynum].push_back(qkappa);
                layered_qshear[laynum].push_back(qshear);
                layered_vph[laynum].push_back(vph / this->velocity_norm);
                layered_vsh[laynum].push_back(vsh / this->velocity_norm);
                layered_eta[laynum].push_back(eta);

                // check whether solid or fluid
                if (idxouter == 0) {
                    if (vsv == 0.0 && vsh == 0.0) {
                        _issolid[laynum] = false;
                    }
                }
            }

            // move to next line
            modelfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            // increment counters
            ++idxinner;
            ++idxouter;
            // }
        }
        // for (auto &idxouter : layered_radii) {
        //    // std::cout << "HELLO\n";
        //    std::cout << idxouter.front() << " " << idxouter.back() << "\n";
        // }
        // std::cout << "Seg test 1\n";
        _numlayers = layered_radii.size();
        vec_layers.reserve(_numlayers + 1);
        vec_layers.push_back(0.0);
        // std::generate(
        //     vec_layers.begin() + 1, vec_layers.end(),
        //     [n = 0, this]() mutable { return layered_radii[n++].back(); });
        // std::cout << "Seg test 2\n";
        _vec_indices =
            std::vector<std::vector<int>>(_numlayers, std::vector<int>(2, 0));
        _vec_indices[0][0] = 0;
        // std::cout << "Seg test 3\n";
        for (int idx = 0; idx < _numlayers; ++idx) {
            vec_layers.push_back(layered_radii[idx].back());
            if (idx != 0) {
                _vec_indices[idx][0] = _vec_indices[idx - 1][1] + 1;
            }
            _vec_indices[idx][1] =
                _vec_indices[idx][0] + layered_radii[idx].size() - 1;
            // std::cout << "Seg test " << idx + 3 << "\n";
        }
        // std::cout << "Size: " << vec_layers.size() << "\n";
        // for (auto &idx : vec_layers) {
        //    std::cout << idx << "\n";
        // }

        for (int idx = 0; idx < _numlayers; ++idx) {
            // iterators to start and end of layer of radius
            auto it1 = layered_radii[idx].begin();
            auto it2 = layered_radii[idx].end();

            // iterators to beginning of this layer for all data
            auto it_rho = layered_rho[idx].begin();
            auto it_vpv = layered_vpv[idx].begin();
            auto it_vsv = layered_vsv[idx].begin();
            auto it_qkappa = layered_qkappa[idx].begin();
            auto it_qshear = layered_qshear[idx].begin();
            auto it_vph = layered_vph[idx].begin();
            auto it_vsh = layered_vsh[idx].begin();
            auto it_eta = layered_eta[idx].begin();

            // pushback
            func_rho.push_back(InterpA(it1, it2, it_rho));
            func_vpv.push_back(InterpA(it1, it2, it_vpv));
            func_vsv.push_back(InterpA(it1, it2, it_vsv));
            func_qkappa.push_back(InterpA(it1, it2, it_qkappa));
            func_qshear.push_back(InterpA(it1, it2, it_qshear));
            func_vph.push_back(InterpA(it1, it2, it_vph));
            func_vsh.push_back(InterpA(it1, it2, it_vsh));
            func_eta.push_back(InterpA(it1, it2, it_eta));
        }

        // for (int idx = 0; idx < N; ++idx) {
        //    double radius, rho, vpv, vsv, qkappa, qshear, vph, vsh, eta;
        //    modelfile >> radius >> rho >> vpv >> vsv >> qkappa >> qshear >>
        //    vph
        //    >>
        //        vsh >> eta;
        //    vec_radius.push_back(radius);
        //    // vec_rho.push_back(rho);
        //    // vec_vpv.push_back(vpv);
        //    // vec_vsv.push_back(vsv);
        //    // vec_qkappa.push_back(qkappa);
        //    // vec_qshear.push_back(qshear);
        //    // vec_vph.push_back(vph);
        //    // vec_vsh.push_back(vsh);
        //    // vec_eta.push_back(eta);

        //    // move to next line
        //    modelfile.ignore(std::numeric_limits<std::streamsize>::max(),
        //    '\n');
        // }

        modelfile.close();
    } else {
        assert("Model not found!");
    }

    // finding layers
    // this->vec_layers = this->findlayers(this->vec_radius);

    // auto vec_indices = this->layerindices(this->vec_radius);
    // for (auto &idx : vec_indices) {
    //    std::cout << idx << "\n";
    // }
    // {
    //    int idxouter = 0;
    //    for (int idxlayers = 0; idxlayers < _numlayers; ++idxlayers) {
    //       std::vector<double> tmp;
    //       while (vec_radius[idxouter] != vec_layers[idxlayers + 1]) {
    //          tmp.push_back(vec_radius[idxouter]);
    //       }
    //    }
    // }
};

};   // namespace EarthModels

#endif