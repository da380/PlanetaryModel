#ifndef PLANETARY_MODEL_CONCEPTS_GUARD_H
#define PLANETARY_MODEL_CONCEPTS_GUARD_H

#include <concepts>
#include <ranges>

namespace PlanetaryModel {

// Concept for a spherical geometry model.
template <typename Model>
concept SphericalGeometryModel = requires(Model model, Model::size_type i,
                                          Model::value_type r) {
  // Member type "size_type" needed that must be integral.
  typename Model::size_type;
  requires std::integral<typename Model::size_type>;

  // Member type "value_type" needed that must be floating point.
  typename Model::value_type;
  requires std::floating_point<typename Model::value_type>;

  // Member function returning the number of layers in the model.
  { model.NumberOfLayers() } -> std::convertible_to<typename Model::size_type>;

  // Member function returning the length normalisation.
  { model.LengthNorm() } -> std::convertible_to<typename Model::value_type>;

  // Member function returning the radius of the ith layer.
  { model.Radius(i) } -> std::convertible_to<typename Model::value_type>;

  // Member function returning the outer radius.
  { model.OuterRadius() } -> std::convertible_to<typename Model::value_type>;

  // Member function returning view to the radii.
  { model.Radii() } -> std::ranges::common_range<>;
};

// Concept for a spherical density model.
template <typename Model>
concept SphericalDensityModel = requires(Model model, Model::size_type i,
                                         Model::value_type r) {
  // Needs all properties of a spherical geometry model.
  requires SphericalGeometryModel<Model>;

  // Member function returning the mass normalisation.
  { model.MassNorm() } -> std::convertible_to<typename Model::value_type>;

  // Member function returning the density normalisation.
  { model.DensityNorm() } -> std::convertible_to<typename Model::value_type>;

  // Member function returning the moment of Inertia normalisation.
  { model.InertiaNorm() } -> std::convertible_to<typename Model::value_type>;

  // Member function returning the time normalisation.
  { model.TimeNorm() } -> std::convertible_to<typename Model::value_type>;

  // Member function returning the velocity normalisation.
  { model.VelocityNorm() } -> std::convertible_to<typename Model::value_type>;

  // Member function returning the acceleration normalisation.
  {
    model.AccelerationNorm()
    } -> std::convertible_to<typename Model::value_type>;

  // Member function returning the force normalisation.
  { model.ForceNorm() } -> std::convertible_to<typename Model::value_type>;

  // Member function returning the stress normalisation.
  { model.StressNorm() } -> std::convertible_to<typename Model::value_type>;

  // Member function returning the gravitational constant.
  {
    model.GravitationalConstant()
    } -> std::convertible_to<typename Model::value_type>;

  // Member function to return density in the ith layer.
  { model.Density(i) } -> std::regular_invocable<typename Model::value_type>;
  { model.Density(i)(r) } -> std::convertible_to<typename Model::value_type>;

  // Member function to return density in ith layer at given radius.
  { model.Density(i, r) } -> std::convertible_to<typename Model::value_type>;
};

template <typename Model>
concept SphericalElasticModel = requires(Model model, Model::size_type i,
                                         Model::value_type r) {
  // Needs all properties of a spherical density model.
  requires SphericalDensityModel<Model>;
};

}  // namespace PlanetaryModel

#endif  // PLANETARY_MODEL_CONCEPTS_GUARD_H
