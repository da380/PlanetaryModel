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

    // define domain:
    double rmax = 10.0;
    double rmin = 0.0;
    int nelem = 10;
    int npoly = 5;
    int matlen = nelem * npoly + 1;
    double dr = (rmax - rmin) / nelem;
    Eigen::MatrixXd matspec = Eigen::MatrixXd::Zero(matlen, matlen);

    Eigen::MatrixXd matsub1 = Eigen::MatrixXd::Zero(npoly + 1, npoly + 1);
    Eigen::MatrixXd matsub2 = Eigen::MatrixXd::Zero(npoly + 1, npoly + 1);
    Eigen::VectorXd vecforce = Eigen::VectorXd::Zero(matlen);
    Eigen::VectorXd vecsol = Eigen::VectorXd::Zero(matlen);
    // find each sub-block for mass matrix
    // generate Gauss grid
    auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);
    std::vector<double> vval(nelem, 1.0);

    auto pleg = LagrangePolynomial(q.Points().begin(), q.Points().end());
    auto func = [pleg](int i, int j, double xval) {
        return pleg(i, xval) * pleg(j, xval);
    };
    auto func1 = [pleg](int i, int j, double xval) {
        return pleg.Derivative(i, xval) * pleg(j, xval);
    };

    for (int i = 0; i < npoly + 1; ++i) {
        for (int j = 0; j < npoly + 1; ++j) {
            auto funcij = [func, i, j](double xval) {
                return func(i, j, xval);
            };
            auto funcij2 = [func1, i, j](double xval) {
                return func1(i, j, xval);
            };
            matsub1(i, j) = q.Integrate(funcij);
            matsub2(i, j) = q.Integrate(funcij2);
        }
    }
    // force vector
    for (int idx = 0; idx < matlen; ++idx) {
        auto funcforce = [pleg, idx, npoly](double x) {
            return pleg(idx % npoly, x);
        };
        vecforce(idx) += q.Integrate(funcforce);
    }
    // std::cout << matsub1 << std::endl;
    // std::cout << matsub2 << std::endl;

    // std::cout << "\n Value of 0th Lagrange polynomial and derivative at -1: "
    //           << std::endl;
    // std::cout << pleg(0, -1) << std::endl;
    // std::cout << pleg.Derivative(0, -1) << std::endl;

    // fill out matspec
    for (int idx = 0; idx < nelem; ++idx) {
        matspec.block(idx * npoly, idx * npoly, npoly + 1, npoly + 1) -=
            matsub2;
    }

    // final element
    matspec(matlen - 1, matlen - 1) += pleg(npoly, 1.0) * pleg(npoly, 1.0);
    Eigen::FullPivLU<Eigen::MatrixXd> solver(matspec);
    vecsol = solver.solve(vecforce);
    auto file = std::ofstream("SpecElem.out");
    for (int i = 0; i < matlen; ++i) {
        file << vecsol(i) << std::endl;
    }
    // return
    // int n = 100;
    // dx = (2.0) / static_cast<double>(n - 1);
    // auto file = std::ofstream("Lagrange.out");
    // for (int i = 0; i < n; i++) {
    //     double x = -1.0 + i * dx;
    //     int j = 1;
    //     file << x << ";" << pleg(0, x) << ";" << pleg.Derivative(0, x)
    //          << std::endl;
    // }
}
