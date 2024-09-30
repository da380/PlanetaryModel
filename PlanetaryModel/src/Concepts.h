#ifndef PLANETARY_MODEL_CONCEPTS_GUARD_H
#define PLANETARY_MODEL_CONCEPTS_GUARD_H

#include <concepts>
#include <ranges>

namespace PlanetaryModel {

template <typename T>
concept HasNormalisationInformation = requires(T t) {
    // Member type "value_type" needed that must be floating point.
    typename T::value_type;
    requires std::floating_point<typename T::value_type>;

    // Member functions for returning base normalisations.
    { t.LengthNorm() } -> std::convertible_to<typename T::value_type>;
    { t.MassNorm() } -> std::convertible_to<typename T::value_type>;
    { t.TimeNorm() } -> std::convertible_to<typename T::value_type>;

    // Member functions for returning derived normalisations.
    { t.DensityNorm() } -> std::convertible_to<typename T::value_type>;
    { t.InertiaNorm() } -> std::convertible_to<typename T::value_type>;
    { t.VelocityNorm() } -> std::convertible_to<typename T::value_type>;
    { t.AccelerationNorm() } -> std::convertible_to<typename T::value_type>;
    { t.ForceNorm() } -> std::convertible_to<typename T::value_type>;
    { t.StressNorm() } -> std::convertible_to<typename T::value_type>;
    {
        t.GravitationalConstant()
    } -> std::convertible_to<typename T::value_type>;
};

// Concept for a spherical geometry model.
template <typename Model>
concept SphericalGeometryModel =
    requires(Model model, Model::size_type i, Model::value_type r) {
        // Member type "size_type" needed that must be integral.
        typename Model::size_type;
        requires std::integral<typename Model::size_type>;

        // Member type "value_type" needed that must be floating point.
        typename Model::value_type;
        requires std::floating_point<typename Model::value_type>;

        // Provides normalisations.
        requires HasNormalisationInformation<Model>;

        // Member function returning the number of layers in the model.
        {
            model.NumberOfLayers()
        } -> std::convertible_to<typename Model::size_type>;

        // Member functions returning the lower and upper radius of the ith
        // layer.
        {
            model.LowerRadius(i)
        } -> std::convertible_to<typename Model::value_type>;
        {
            model.UpperRadius(i)
        } -> std::convertible_to<typename Model::value_type>;

        // Member function returning the outer radius.
        {
            model.OuterRadius()
        } -> std::convertible_to<typename Model::value_type>;
    };

// Concept for a spherical density model.
template <typename Model>
concept SphericalDensityModel =
    requires(Model model, Model::size_type i, Model::value_type r) {
        // Needs all properties of a spherical geometry model.
        requires SphericalGeometryModel<Model>;

        // Member function to return density in the ith layer.
        {
            model.Density(i)
        } -> std::regular_invocable<typename Model::value_type>;
        {
            model.Density(i)(r)
        } -> std::convertible_to<typename Model::value_type>;
    };

template <typename Model>
concept SphericalElasticModel =
    requires(Model model, Model::size_type i, Model::value_type r) {
        // Needs all properties of a spherical density model.
        requires SphericalDensityModel<Model>;

        // Member function to determine is model is isotropic.
        { model.IsIsotropic() } -> std::same_as<bool>;

        // Member functions to determine whether a layer is solid or fluid.
        { model.IsSolid(i) } -> std::same_as<bool>;
        { model.IsFluid(i) } -> std::same_as<bool>;

        // Member functions to return TI elastic modulii.
        { model.A(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.A(i)(r) } -> std::convertible_to<typename Model::value_type>;

        { model.C(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.C(i)(r) } -> std::convertible_to<typename Model::value_type>;

        { model.F(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.F(i)(r) } -> std::convertible_to<typename Model::value_type>;

        { model.L(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.L(i)(r) } -> std::convertible_to<typename Model::value_type>;

        { model.N(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.N(i)(r) } -> std::convertible_to<typename Model::value_type>;

        // Member functions to return isotropic elastic modulii.
        {
            model.Kappa(i)
        } -> std::regular_invocable<typename Model::value_type>;
        {
            model.Kappa(i)(r)
        } -> std::convertible_to<typename Model::value_type>;

        { model.Mu(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.Mu(i)(r) } -> std::convertible_to<typename Model::value_type>;

        // Member functions to return TI elastic wave speeds.
        { model.VPV(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.VPV(i)(r) } -> std::convertible_to<typename Model::value_type>;

        { model.VPH(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.VPH(i)(r) } -> std::convertible_to<typename Model::value_type>;

        { model.VSV(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.VSV(i)(r) } -> std::convertible_to<typename Model::value_type>;

        { model.VSH(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.VSH(i)(r) } -> std::convertible_to<typename Model::value_type>;

        // Member functions to return isotropic elastic wave speeds.
        { model.VP(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.VP(i)(r) } -> std::convertible_to<typename Model::value_type>;

        { model.VS(i) } -> std::regular_invocable<typename Model::value_type>;
        { model.VS(i)(r) } -> std::convertible_to<typename Model::value_type>;
    };

// Concept for a spherical density model.
template <typename Model>
concept AsphericalDensityModel =
    requires(Model model, Model::size_type i, Model::value_type r,
             Model::value_type theta, Model::value_type phi) {
        // Needs all properties of a spherical geometry model.
        requires SphericalDensityModel<Model>;

        // Member function to return density in the ith layer.
        {
            model.DensityPerturbation(i)
        } -> std::regular_invocable<typename Model::value_type>;
        {
            model.DensityPerturbation(i)(r, theta, phi)
        } -> std::convertible_to<typename Model::value_type>;
    };

// Concept for a simple deviation model.
template <typename Model>
concept SingleParameterDeviationModel =
    requires(Model model, int i, double depth, double theta, double phi) {
        // Member function to return density in the ith layer.
        { model.GetValueAt(depth, theta, phi) } -> std::convertible_to<double>;
    };
}   // namespace PlanetaryModel

#endif   // PLANETARY_MODEL_CONCEPTS_GUARD_H
