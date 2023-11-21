#include <GaussQuad/All>
#include <PlanetaryModel/All>
// #include <boost/units/systems/si/current.hpp>
// #include <boost/units/systems/si/electric_potential.hpp>
// #include <boost/units/systems/si/energy.hpp>
// #include <boost/units/systems/si/force.hpp>
// #include <boost/units/systems/si/io.hpp>
// #include <boost/units/systems/si/length.hpp>
// #include <boost/units/systems/si/resistance.hpp>

#include <complex>
#include <concepts>
#include <iostream>

#include "PREM.h"
#include <boost/mpl/list.hpp>

#include <boost/typeof/std/complex.hpp>

#include <boost/units/io.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/quantity.hpp>

#include "test_system.hpp"
int
main() {
    using Float = double;

    using namespace GaussQuad;
    using namespace boost::units;
    using namespace boost::units::test;
    // Set the output precision
    using std::cout;
    using std::endl;
    cout.setf(std::ios_base::scientific);
    cout.setf(std::ios_base::showpos);
    cout.precision(16);

    // Set the floating point precision
    using Float = double;

    // Build the quadrature
    int n = 5;
    auto q = GaussLobattoLegendreQuadrature1D<Float>(n);

    // write out the points and weights
    for (int i = 0; i < n; i++) {
        std::cout << q.X(i) << " " << q.W(i) << std::endl;
    }

    // define a simple function to integrate
    auto fun = [](Float x) { return x * x; };

    // set the exact value for the integral
    Float exact = Float(2.0) / Float(3.0);

    cout << "Numerical value = " << q.Integrate(fun)
         << ", exact value = " << exact << endl;

    // assertions of concept for model
    static_assert(
        PlanetaryModel::HasNormalisationInformation<EarthConstants<double>>);
    static_assert(PlanetaryModel::SphericalGeometryModel<PREM<double, int>>);
    static_assert(PlanetaryModel::SphericalDensityModel<PREM<double, int>>);
    static_assert(PlanetaryModel::SphericalElasticModel<PREM<double, int>>);

    // finding mass and moment of inertia of a spherical model
    auto myprem = PREM();
    std::cout << "Oh Hello: " << myprem.A(0)(0.0) << std::endl;
    double PMass = 0;
    double PINIT = 0;
    int numlayers = myprem.NumberOfLayers();
    auto nmap = [](Float xi, Float x1, Float x2) {
        return 0.5 * ((x1 + x2) + (x2 - x1) * xi);
    };
    auto dval = [&myprem, nmap](Float x, int i) {
        return myprem.Density(i)(
                   nmap(x, myprem.LowerRadius(i), myprem.UpperRadius(i)) /
                   myprem.OuterRadius()) *
               nmap(x, myprem.LowerRadius(i), myprem.UpperRadius(i)) *
               nmap(x, myprem.LowerRadius(i), myprem.UpperRadius(i));
    };
    auto IVAL = [&myprem, nmap](Float x, int i) {
        return myprem.Density(i)(
                   nmap(x, myprem.LowerRadius(i), myprem.UpperRadius(i)) /
                   myprem.OuterRadius()) *
               std::pow(nmap(x, myprem.LowerRadius(i), myprem.UpperRadius(i)),
                        4);
    };

    std::cout << "Hello: " << nmap(-1, 0, 10) << " " << nmap(1, 0, 10)
              << std::endl;
    std::cout << "Hello: " << dval(-1.0, 0) << " " << dval(1.0, 0) << std::endl;
    std::cout << myprem.Density(0)(
                     nmap(-1, myprem.LowerRadius(0), myprem.UpperRadius(0)))
              << std::endl;

    // loop over all layers
    for (int idx = 0; idx < myprem.NumberOfLayers(); ++idx) {
        auto lval = [dval, idx](Float x) { return dval(x, idx); };
        PMass += 1000.0 * 4.0 * 3.1415926535 * 0.5 *
                 (myprem.UpperRadius(idx) - myprem.LowerRadius(idx)) *
                 q.Integrate(lval);
        auto lval2 = [IVAL, idx](Float x) { return IVAL(x, idx); };
        PINIT += 1000.0 * 8.0 * 3.1415926535 / 3 * 0.5 *
                 (myprem.UpperRadius(idx) - myprem.LowerRadius(idx)) *
                 q.Integrate(lval2);
    }
    std::cout << "The mass of the Earth in PREM is: " << PMass << " kg"
              << std::endl;
    std::cout << "The moment of inertia of the Earth in PREM is: " << PINIT
              << std::endl;

    {
        //[quantity_snippet_1
        quantity<length> L = 2.0 * meters;   // quantity of length
        quantity<energy> E =
            kilograms * pow<2>(L / seconds);   // quantity of energy
        //]

        std::cout
            << "L                                 = " << L << std::endl
            << "L+L                               = " << L + L << std::endl
            << "L-L                               = " << L - L << std::endl
            << "L*L                               = " << L * L << std::endl
            << "L/L                               = " << L / L << std::endl
            << "L*meter                           = " << L * meter << std::endl
            << "kilograms*(L/seconds)*(L/seconds) = "
            << kilograms * (L / seconds) * (L / seconds) << std::endl
            << "kilograms*(L/seconds)^2           = "
            << kilograms * pow<2>(L / seconds) << std::endl
            << "L^3                               = " << pow<3>(L) << std::endl
            << "L^(3/2)                           = "
            << pow<static_rational<3, 2>>(L) << std::endl
            << "2vL                               = " << root<2>(L) << std::endl
            << "(3/2)vL                           = "
            << root<static_rational<3, 2>>(L) << std::endl
            << std::endl;
    }
}
