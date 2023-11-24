#include <Eigen/Core>
#include <Eigen/Dense>
#include <GaussQuad/All>
#include <Interpolation/All>
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
#include <fstream>
#include <iostream>

#include "EARTHMODELS.h"
#include <boost/mpl/list.hpp>

#include <boost/typeof/std/complex.hpp>

#include <boost/units/io.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/quantity.hpp>

#include "test_system.hpp"
int
main() {
    using namespace Interpolation;

    // // Set the nodes
    // int N = 3;
    // double x1 = 0.0;
    // double x2 = 1.0;
    // double dx = (x2 - x1) / static_cast<double>(N - 1);
    // std::vector<double> X(N);
    // for (int i = 0; i < N; i++) {
    //     X[i] = x1 + i * dx;
    // }

    // // Set the Lagrange polynomial
    // auto p = LagrangePolynomial(X.begin(), X.end());

    // // Set the values for plotting
    // int n = 100;
    // dx = (x2 - x1) / static_cast<double>(n - 1);

    // auto file = std::ofstream("Lagrange.out");
    // for (int i = 0; i < n; i++) {
    //     double x = x1 + i * dx;
    //     int j = 1;
    //     file << x << " " << p(j, x) << " " << p.Derivative(j, x) <<
    //     std::endl;
    // }

    // function:
    int N = 5;
    std::vector<double> x(N), y(N);
    // int N = 3;
    double x1 = 0.0;
    double x2 = 1.0;
    double dx = (x2 - x1) / static_cast<double>(N - 1);
    for (int idx = 0; idx < N; ++idx) {
        x[idx] = dx * idx;
    };
    // std::generate(y.begin(), y.end(),
    //               [nidx, x]() { return x[nidx] * x[nidx]; });
    std::transform(x.begin(), x.end(), y.begin(),
                   [](double x) { return x * x; });
    for (int idx = 0; idx < N; ++idx) {
        std::cout << x[idx] << "  " << y[idx] << std::endl;
    }
    auto plag = Lagrange(x.begin(), x.end(), y.begin());

    // return
    int n = 100;
    dx = (x2 - x1) / static_cast<double>(n - 1);
    auto file = std::ofstream("Lagrange.out");
    for (int i = 0; i < n; i++) {
        double x = x1 + i * dx;
        int j = 1;
        file << x << " " << plag(x) << " " << plag.Derivative(x) << std::endl;
    }

    // define domain:
    double rmax = 10.0;
    double rmin = 0.0;
    int nelem = 10;
    int npoly = 5;
    int matlen = nelem * npoly + 1;
    double dr = (rmax - rmin) / nelem;
    Eigen::MatrixXd matspec = Eigen::MatrixXd::Zero(matlen, matlen);
    Eigen::MatrixXd matsub = Eigen::MatrixXd::Zero(npoly + 1, npoly + 1);

    // find each sub-block for mass matrix
    // generate Gauss grid
    auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);
    std::vector<double> vval(nelem, 1.0);

    auto pleg = LagrangePolynomial(q.Points().begin(), q.Points().end());
    auto func = [pleg](int i, int j, double xval) {
        return pleg(i, xval) * pleg(j, xval);
    };

    for (int i = 0; i < npoly + 1; ++i) {
        for (int j = 0; j < npoly + 1; ++j) {
            auto funcij = [func, i, j](double xval) {
                return func(i, j, xval);
            };
            matsub(i, j) = q.Integrate(funcij);
        }
    }
    std::cout << matsub << std::endl;
}
