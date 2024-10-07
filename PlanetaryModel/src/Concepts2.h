#ifndef TOMOGRAPHY_MODEL_CONCEPTS_GUARD_H
#define TOMOGRAPHY_MODEL_CONCEPTS_GUARD_H

#include <concepts>
#include <ranges>
#include <vector>

namespace PlanetaryModel {

// Concept for a spherical geometry model.
template <typename Model>
concept TomographyModel = requires(Model model, double d, double lo, double la,
                                   bool outofBound = false) {
    // value of tomography model
    { model.GetValueAt(d, lo, la, outofBound) } -> std::convertible_to<double>;

    // depths
    { model.GetDepths() } -> std::convertible_to<std::vector<double>>;

    // Latitudes
    { model.GetLatitudes() } -> std::convertible_to<std::vector<double>>;

    // Longitudes
    { model.GetLongitudes() } -> std::convertible_to<std::vector<double>>;

    // Values
    { model.GetValues() } -> std::convertible_to<std::vector<double>>;
};

// zero change tomography model
class TomographyZeroModel {
  public:
    double GetValueAt(double d, double lo, double la, bool outofBound = false) {
        return 0.0;
    };
    auto GetDepths() { return vec_depths; };
    auto GetLatitudes() { return vec_latitudes; };
    auto GetLongitudes() { return vec_longitudes; };
    auto GetValues() { return vec_values; };

  private:
    std::vector<double> vec_depths{0, 6371.0}, vec_latitudes{-90.0, 90.0},
        vec_longitudes{0.0, 360.0}, vec_values{0};
};

}   // namespace PlanetaryModel

#endif   // PLANETARY_MODEL_CONCEPTS_GUARD_H
