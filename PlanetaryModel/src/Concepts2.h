#ifndef TOMOGRAPHY_MODEL_CONCEPTS_GUARD_H
#define TOMOGRAPHY_MODEL_CONCEPTS_GUARD_H

#include <concepts>
#include <ranges>

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

}   // namespace PlanetaryModel

#endif   // PLANETARY_MODEL_CONCEPTS_GUARD_H
