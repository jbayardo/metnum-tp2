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
#include "Matrix.h"

template <typename T>
using min_queue = std::priority_queue<T, std::vector<T>, std::greater<T>>;
// Matriz de datos, vector, número de línea.
typedef std::function<double(const Matrix &, const std::vector<double> &, int)> DistanceF;
typedef std::function<double(const std::vector<double> &)> Norm;
typedef std::function<bool(const Matrix &, const std::vector<double> &, const double &, unsigned int)> ConditionF;
typedef std::pair<double, std::vector<double>> EigenPair;

std::vector<double> operator*(const Matrix &m, const std::vector<double> &n);
EigenPair powerIteration(const Matrix &A, std::vector<double> eigenVector, const Norm &norm, const ConditionF &condition);
Matrix deflation(const Matrix &A, const EigenPair &eigen);
std::list<EigenPair> decompose(const Matrix &A, int k, const Norm &norm, const ConditionF &condition);

const DistanceF L2 = DistanceF([](const Matrix &A, const std::vector<double> &v, int i) -> double {
    if (v.size() < A.columns()) {
        throw new std::out_of_range("Invalid size for vector");
    }

    double output = 0.0;

    for (int j = 0; j < A.columns(); ++j) {
        output += std::pow(A(i, j) - v[i], 2.0);
    }

    return std::sqrt(output);
});

const Norm N2 = Norm([](const std::vector<double> &v) -> double {
    return std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
});

#endif //METNUM_TP2_LINEARALGEBRA_H
