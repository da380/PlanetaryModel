#ifndef TOMOGRAPHY_MODEL_CONCEPTS_GUARD_H
#define TOMOGRAPHY_MODEL_CONCEPTS_GUARD_H

#include <concepts>
#include <ranges>
#include <vector>

namespace PlanetaryModel {

// simplified concepts
template <typename T>
concept HasBasicNormalisationInformation = requires(T t) {
    // Member functions for returning base normalisations.
    { t.LengthNorm() } -> std::convertible_to<double>;
    { t.MassNorm() } -> std::convertible_to<double>;
    { t.TimeNorm() } -> std::convertible_to<double>;
};

// Concept for a spherical geometry model.
template <typename Model, typename INT, typename FLOAT>
concept BasicSphericalGeometryModel = requires(Model model, INT i, FLOAT r) {
    // Member type "size_type" needed that must be integral.
    // typename Model::size_type;
    requires std::integral<INT>;

    // Member type "value_type" needed that must be floating point.
    // typename Model::value_type;
    requires std::floating_point<FLOAT>;

    // Provides normalisations.
    requires HasBasicNormalisationInformation<Model>;

    // Member function returning the number of layers in the model.
    { model.NumberOfLayers() } -> std::convertible_to<INT>;

    // Member functions returning the lower and upper radius of the ith
    // layer.
    { model.LowerRadius(i) } -> std::convertible_to<FLOAT>;
    { model.UpperRadius(i) } -> std::convertible_to<FLOAT>;

    // Member function returning the outer radius.
    { model.OuterRadius() } -> std::convertible_to<FLOAT>;
};

// Concept for a spherical density model.
template <typename Model, typename INT, typename FLOAT>
concept BasicSphericalDensityModel = requires(Model model, INT i, FLOAT r) {
    // Needs all properties of a spherical geometry model.
    requires BasicSphericalGeometryModel<Model, INT, FLOAT>;

    // Member function to return density in the ith layer.
    { model.Density(i) } -> std::regular_invocable<FLOAT>;
    { model.Density(i)(r) } -> std::convertible_to<FLOAT>;
};

// Concept for a spherical density model.
template <typename Model, typename INT, typename FLOAT>
concept BasicAsphericalDensityModel =
    requires(Model model, INT i, FLOAT r, FLOAT theta, FLOAT phi) {
        // Needs all properties of a spherical geometry model.
        requires BasicSphericalDensityModel<Model, INT, FLOAT>;

        // Member function to return density in the ith layer.
        { model.Mapping(i) } -> std::regular_invocable<FLOAT, FLOAT, FLOAT>;
        { model.Mapping(i)(r, theta, phi) } -> std::convertible_to<FLOAT>;
    };

// concept for a mapping class
template <typename Model>
concept MappingClass =
    requires(Model model, int i, double r, double theta, double phi) {
        // Member function to return density in the ith layer.
        { model.Mapping(i) } -> std::regular_invocable<double, double, double>;
        { model.Mapping(i)(r, theta, phi) } -> std::convertible_to<double>;
    };
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
    TomographyZeroModel() {};

    const double GetValueAt(double d, double lo, double la,
                            bool outofBound = false) const {
        return 0.0;
    };
    const auto GetDepths() const { return vec_depths; };
    const auto GetLatitudes() const { return vec_latitudes; };
    const auto GetLongitudes() const { return vec_longitudes; };
    const auto GetValues() const { return vec_values; };

  private:
    const std::vector<double> vec_depths{0, 6371.0}, vec_latitudes{-90.0, 90.0},
        vec_longitudes{0.0, 360.0}, vec_values{0};
};

}   // namespace PlanetaryModel

#endif   // PLANETARY_MODEL_CONCEPTS_GUARD_H
