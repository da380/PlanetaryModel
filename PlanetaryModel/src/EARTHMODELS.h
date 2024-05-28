#ifndef EARTHMODELS_GUARD_H
#define EARTHMODELS_GUARD_H

#include <PlanetaryModel/All>
#include <cmath>
#include <concepts>
#include <functional>
#include <iostream>
#include <ranges>
#include <vector>

#include "Interpolation/All"

namespace EarthModels {
template <typename FLOAT> class EarthConstants {
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
    std::vector<FLOAT> vec_radii{0.0,       1221500.0, 3480000.0, 3630000.0,
                                 5600000.0, 5701000.0, 5771000.0, 5971000.0,
                                 6151000.0, 6291000.0, 6346600.0, 6356000.0,
                                 6368000.0, 6371000.0};

    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_density{
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
    PERTPREM(){};

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
    HOMOBOUND0(){};

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
    HOMOBOUND1(){};

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
    HOMOBOUND2(){};

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
    HOMOBOUND3(){};

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
    HOMOBOUND4(){};

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
    HOMOBOUND5(){};

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
    HOMOBOUND6(){};

    // Density
    Interpolation::Polynomial1D<FLOAT> DensityPerturbation(INTEGRAL i) {
        return vec_pert_density[i];
    };
    // std::function<FLOAT(FLOAT, FLOAT, FLOAT)> RadialMap() const {
    //     return RadialMap();
    // };
    FLOAT RadialMap(FLOAT r, FLOAT theta, FLOAT phi) const {
        return 0.01 * r * std::sin(theta);
    };
    FLOAT MaxRadius() const { return this->OuterRadius() * 1.01; };

  private:
    std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pert_density{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

};   // namespace EarthModels

#endif