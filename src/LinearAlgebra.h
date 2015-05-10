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
#include "Matrix.h"

template <typename T>
using min_queue = std::priority_queue<T, std::vector<T>, std::greater<T>>;
// Matriz de datos, vector, número de línea.
typedef std::function<double(const Matrix &, int, const Matrix &, int)> DistanceF;
typedef std::function<double(const std::vector<double> &)> Norm;
typedef std::function<bool(const Matrix &, const std::vector<double> &, const double &, unsigned int)> ConditionF;
typedef std::pair<double, std::vector<double>> EigenPair;

std::vector<double> operator*(const Matrix &m, const std::vector<double> &n);
EigenPair powerIteration(const Matrix &A, std::vector<double> eigenVector, const Norm &norm, const ConditionF &condition);
void deflation(const Matrix &A, const EigenPair &eigen);
std::list<EigenPair> decompose(const Matrix &A, int k, const Norm &norm, const ConditionF &condition);
void dimensionReduction(const Matrix& src, Matrix& dst, const std::list<EigenPair>& l);

const DistanceF L2 = DistanceF([](const Matrix &A, int i0, const Matrix &B, int i1) -> double {
    if (B.columns() != A.columns()) {
        throw new std::out_of_range("Invalid size for vector");
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

#endif //METNUM_TP2_LINEARALGEBRA_H
