//
// Created by Julian Bayardo on 5/8/15.
//

#ifndef METNUM_TP2_LINEARALGEBRA_H
#define METNUM_TP2_LINEARALGEBRA_H

#include <vector>
#include <cmath>
#include <queue>
#include <functional>
#include <list>
#include <numeric>
#include <algorithm>
#include <sstream>
#include "Matrix.h"

template <typename T>
using min_queue = std::priority_queue<T, std::vector<T>, std::greater<T>>;
// Matriz de datos, vector, número de línea.
typedef std::function<double(const Matrix &, int, const Matrix &, int)> DistanceF;
typedef std::function<double(const std::vector<double> &)> Norm;
typedef std::function<bool(const Matrix &, const std::vector<double> &, const double &, unsigned int)> ConditionF;
typedef std::pair<double, std::vector<double>> EigenPair;

std::vector<double> operator*(const Matrix &m, const std::vector<double> &n);
std::vector<double> operator*(const double &m, const std::vector<double> &n);
EigenPair powerIteration(const Matrix &, std::vector<double>, const Norm &, const ConditionF &);
void deflation(const Matrix &A, const EigenPair &eigen);
std::list<EigenPair> decompose(Matrix, int, const Norm &, const ConditionF &);
void dimensionReduction(const Matrix& src, Matrix& dst, const std::list<EigenPair>& l);

const DistanceF L2 = DistanceF([](const Matrix &A, int i0, const Matrix &B, int i1) -> double {
    if (B.columns() != A.columns()) {
        std::stringstream fmt;
        fmt << "Tamaño de matriz A es " << A.columns() << ", mientras que B es " << B.columns();
        throw new std::out_of_range(fmt.str());
    }

    double output = 0.0;

    for (int j = 0; j < A.columns(); ++j) {
        output += std::pow(A(i0, j) - B(i1, j), 2.0);
    }

    return std::sqrt(output);
});

const Norm N2 = Norm([](const std::vector<double> &v) -> double {
    return std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
});

auto CIterations = [](unsigned int limit) -> ConditionF {
    return ConditionF([limit](const Matrix &A, const std::vector<double> &v, const double &c, unsigned int its) -> bool {
        return its > limit;
    });
};

/*auto CPrecision = [](double epsilon) -> ConditionF {
    return ConditionF([epsilon](const Matrix &A, const std::vector<double> &v, const double &c, unsigned int its) -> bool {
        return A * v - c * v < epsilon;
    });
};*/

#endif //METNUM_TP2_LINEARALGEBRA_H
