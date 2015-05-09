//
// Created by Julian Bayardo on 4/25/15.
//

#ifndef TP2_MATRIX_H
#define TP2_MATRIX_H

#include <iostream>
#include <bitset>

// Este magic number nos dice cuándo convertir automáticamente una matriz banda en una matriz normal.
#define MAGIC_NUMBER 562154

/*
* Matriz Banda.
*/
class Matrix {
    friend std::ostream &operator<<(std::ostream &, const Matrix &);
public:
    Matrix(const Matrix &m);

    template <int K>
    Matrix(const Matrix &m, const std::bitset<K> &filter);

    Matrix(int N, int M, int lband = MAGIC_NUMBER, int uband = MAGIC_NUMBER);

    int inline rows() const {
        return this->N;
    }

    int inline columns() const {
        return this->M;
    }

    int inline upper_bandwidth() const {
        return this->uband;
    }

    int inline lower_bandwidth() const {
        return this->lband;
    }

    double & operator()(const int &i, const int &j);
    const double & operator()(const int &i, const int &j) const;
    Matrix & operator=(const Matrix &m);
    bool operator==(const Matrix &m) const;
    bool operator!=(const Matrix &m) const;
    Matrix & operator+=(const Matrix &m);
    Matrix & operator*=(const double &c);

    ~Matrix();
private:
    // Matrix
    int N;
    int M;
    int uband;
    int lband;
    double **matrix;
};

std::ostream & operator<<(std::ostream &os, const Matrix &m);
Matrix operator+(const Matrix &m, const Matrix &n);
Matrix operator*(const Matrix &m, const double &c);
Matrix operator*(const Matrix &m, const Matrix &n);

#endif //TP2_MATRIX_H
