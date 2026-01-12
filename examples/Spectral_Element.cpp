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

#include <algorithm>
#include <complex>
#include <concepts>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

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
    // for (int idx = 0; idx < N; ++idx) {
    //     std::cout << x[idx] << "  " << y[idx] << std::endl;
    // }
    auto plag = Lagrange(x.begin(), x.end(), y.begin());

    // define domain:
    double rmax;
    std::cout << "Maximum value of x: \n";
    std::cin >> rmax;
    int nelem;
    std::cout << "Number of elements: \n";
    std::cin >> nelem;
    bool relem;
    std::cout << "Random nodes: \n";
    std::cin >> relem;
    double rmin = 0.0;

    int npoly = 8;
    int matlen = nelem * npoly + 1;
    double dr = (rmax - rmin) / static_cast<double>(nelem);

    // Random vector of beginning and end of elements
    std::vector<double> vec_elem(nelem + 1);
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(rmin, rmax);
    auto gen = [&distr, &eng]() { return distr(eng); };
    std::generate(vec_elem.begin(), vec_elem.end(), gen);

    vec_elem[0] = rmin;
    vec_elem[nelem] = rmax;
    for (int idx = 1; idx < nelem; ++idx) {
        // std::cout << static_cast<double>(idx) / static_cast<double>(nelem)
        //           << std::endl;
        if (relem) {
            vec_elem[idx] =
                rmax * static_cast<double>(idx) / static_cast<double>(nelem);
        }
    }
    // vec_elem[1] = 1.8;
    std::sort(vec_elem.begin(), vec_elem.end());
    // vec_elem.insert(vec_elem.begin(), rmin);
    // vec_elem.push_back(rmax);
    // for (auto i : vec_elem) {
    //     std::cout << vec_elem[i] << std::endl;
    // }
    for (auto i = 0; i < nelem + 1; ++i) {
        // std::cout << vec_elem[i] << std::endl;
    }
    // Getting a random double value
    // double random_double = unif(re);

    // finite element matrices
    Eigen::MatrixXd matspec = Eigen::MatrixXd::Zero(matlen, matlen);
    Eigen::MatrixXd matsub1 = Eigen::MatrixXd::Zero(npoly + 1, npoly + 1);
    Eigen::MatrixXd matsub2 = Eigen::MatrixXd::Zero(npoly + 1, npoly + 1);
    Eigen::VectorXd vecforce = Eigen::VectorXd::Zero(matlen);
    Eigen::VectorXd vecsol = Eigen::VectorXd::Zero(matlen);

    // find each sub-block for mass matrix
    // generate Gauss grid
    auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);
    std::vector<double> vval(nelem, 1.0);

    // vector of Lagrange polynomials:
    // std::vector<Interpolation::LagrangePolynomial> veclag;

    auto pleg = LagrangePolynomial(q.Points().begin(), q.Points().end());
    auto posforce = [](double x) { return std::cos(x); };
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
    for (int idxelem = 0; idxelem < nelem; ++idxelem) {
        for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
            auto funcforce = [posforce, pleg, idxelem, idxpoly, npoly,
                              vec_elem](double x) {
                // return pleg(idxpoly, x) *
                //        posforce(
                //            ((vec_elem[idxelem + 1] - vec_elem[idxelem]) * x +
                //             (vec_elem[idxelem + 1] + vec_elem[idxelem])) *
                //            0.5);
                return pleg(idxpoly, x) *
                       posforce(
                           ((vec_elem[idxelem + 1] - vec_elem[idxelem]) * x +
                            (vec_elem[idxelem + 1] + vec_elem[idxelem])) *
                           0.5);
            };
            // std::cout << idxelem * npoly + idxpoly << std::endl;
            vecforce(idxelem * npoly + idxpoly) +=
                0.5 * (vec_elem[idxelem + 1] - vec_elem[idxelem]) *
                q.Integrate(funcforce);
            // std::cout << q.Integrate(funcforce) << std::endl;
            // vecforce(idxelem * npoly + idxpoly) += q.Integrate(funcforce);
        }
    }

    // for (int idx = 0; idx < matlen; ++idx) {
    //     auto funcforce = [posforce, pleg, idx, npoly, vec_elem](double x) {
    //         return pleg(idx % npoly, x) *
    //                posforce(((vec_elem[idx + 1] - vec_elem[idx]) * x +
    //                          (vec_elem[idx + 1] + vec_elem[idx])) *
    //                         0.5);
    //     };
    //     vecforce(idx) +=
    //         0.5 * (vec_elem[idx + 1] - vec_elem[idx]) *
    //         q.Integrate(funcforce);
    // }
    // std::cout << matsub1 << std::endl;
    // std::cout << matsub2 << std::endl;

    // std::cout << "\n Value of 0th Lagrange polynomial and derivative at -1: "
    //           << std::endl;
    // std::cout << pleg(0, -1) << std::endl;
    // std::cout << pleg.Derivative(0, -1) << std::endl;

    // for (auto i : vec_elem) {
    //     std::cout << vec_elem[i] << std::endl;
    // }
    // fill out matspec
    for (int idx = 0; idx < nelem; ++idx) {
        matspec.block(idx * npoly, idx * npoly, npoly + 1, npoly + 1) -=
            matsub2;
        // std::cout << 0.5 * (vec_elem[idx + 1] - vec_elem[idx]) << std::endl;
        // matspec.block(idx * npoly, idx * npoly, npoly + 1, npoly + 1) -=
        //     matsub2;
    }

    // final element
    // matspec(matlen - 1, matlen - 1) += 0.5 *
    //                                    (vec_elem[nelem] - vec_elem[nelem -
    //                                    1]) * pleg(npoly, 1.0) *
    //                                    pleg(npoly, 1.0);
    matspec(matlen - 1, matlen - 1) += pleg(npoly, 1.0) * pleg(npoly, 1.0);
    // std::cout << matspec << std::endl;
    // std::cout << vecforce << std::endl;
    Eigen::FullPivLU<Eigen::MatrixXd> solver(matspec);
    vecsol = solver.solve(vecforce);
    auto file = std::ofstream("SpecElem.out");
    for (int i = 0; i < matlen; ++i) {
        file << vecsol(i) << std::endl;
    }

    auto file2 = std::ofstream("SpecOutput.out");
    for (int i = 0; i < nelem + 1; ++i) {
        // file2 << vec_elem[i] << ";" << vec_elem[i] << ";" << vecsol(i *
        // npoly)
        //   << std::endl;
        // file2 << std::setprecision(8) << vec_elem[i] << ";"
        //       << std::pow(vec_elem[i], 2) << ";" << vecsol(i * npoly)
        //       << std::endl;
        file2 << std::setprecision(8) << vec_elem[i] << ";"
              << std::sin(vec_elem[i]) << ";" << vecsol(i * npoly) << std::endl;
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
