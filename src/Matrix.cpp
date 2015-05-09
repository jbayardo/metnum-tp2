//
// Created by Julian Bayardo on 4/25/15.
//

#include "Matrix.h"
#include <stdexcept>

const double zero = 0;

// TODO: revisar las bandas y como se settean en los constructores
Matrix::Matrix(const Matrix &m) : N(m.rows()), M(m.columns()), uband(m.upper_bandwidth()), lband(m.lower_bandwidth()) {
    int bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;
    this->matrix = new double*[this->rows()];

    for (int i = 0; i < this->rows(); ++i) {
        this->matrix[i] = new double[bound];

        for (int j = 0; j < bound; ++j) {
            this->matrix[i][j] = m.matrix[i][j];
        }
    }
}

template <int K>
Matrix::Matrix(const Matrix &m, const std::bitset<K> &filter)
        : N((int)filter.count()), M(m.columns()), uband(m.upper_bandwidth()), lband(std::min(m.lower_bandwidth(), N)) {
    if (K != this->rows()) {
        throw new std::out_of_range("Invalid filter for bitset");
    }

    int last = 0;
    int bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;
    this->matrix = new double*[this->rows()];

    for (int i = 0; i < m.rows(); ++i) {
        if (filter.test((std::size_t)i)) {
            this->matrix[last] = new double[bound];

            for (int j = 0; j < bound; ++j) {
                (*this)(last, j) = m(i, j);
            }

            last++;
        }
    }
}

Matrix::Matrix(int N, int M, int lband, int uband)
        : N(N), M(M), uband(uband), lband(lband) {
    if (this->rows() == 0 || this->columns() == 0) {
        throw new std::out_of_range("Invalid matrix dimension");
    }

    if (lband > N) {
        this->lband = N-1;
    }

    if (uband > M) {
        this->uband = M-1;
    }

    int bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;
    this->matrix = new double*[this->rows()];

    for (int i = 0; i < this->rows(); ++i) {
        this->matrix[i] = new double[bound];

        for (int j = 0; j < bound; ++j) {
            this->matrix[i][j] = 0;
        }
    }
}

double &Matrix::operator()(const int &i, const int &j) {
    if (i >= this->rows() || j >= this->columns()) {
        throw new std::out_of_range("Index access out of range");
    }

    if (i <= j + this->lower_bandwidth() && j <= i + this->upper_bandwidth()) {
        return matrix[i][j - i + this->lower_bandwidth()];
    } else {
        throw new std::out_of_range("Out of modifiable range");
    }
}

const double &Matrix::operator()(const int &i, const int &j) const {
    if (i >= this->rows() || j >= this->columns()) {
        throw new std::out_of_range("Index access out of range");
    }

    if (i > j + this->lower_bandwidth()) {
        return zero;
    } else if (j > i + this->upper_bandwidth()) {
        return zero;
    } else {
        return matrix[i][j - i + this->lower_bandwidth()];
    }
}

Matrix &Matrix::operator=(const Matrix &m) {
    if (*this != m) {
        // Limpiar memoria.
        for (int i = 0; i < this->rows(); ++i) {
            delete[] this->matrix[i];
        }

        delete[] this->matrix;

        // Poner información de la representación interna
        this->N = m.rows();
        this->M = m.columns();
        this->lband = m.lower_bandwidth();
        this->uband = m.upper_bandwidth();

        // Crear matriz nueva
        int bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;
        this->matrix = new double*[m.rows()];

        for (int i = 0; i < this->rows(); ++i) {
            this->matrix[i] = new double[bound];

            for (int j = 0; j < bound; ++j) {
                // Copiar los valores de la matriz
                this->matrix[i][j] = m.matrix[i][j];
            }
        }
    }

    return *this;
}

bool Matrix::operator==(const Matrix &m) const {
    if (this->rows() != m.rows() || this->columns() != m.columns()) {
        return false;
    } else {
        for (int i = 0; i < this->rows(); i++) {
            for (int j = 0; j < this->columns(); j++) {
                if ((*this)(i, j) != m(i, j)) {
                    return false;
                }
            }
        }

        return true;
    }
}

bool Matrix::operator!=(const Matrix &m) const {
    return !(*this == m);
}

Matrix &Matrix::operator+=(const Matrix &m) {
    if (this->rows() == m.rows() && this->columns() == m.columns()) {
        // Si podemos sumar

        if (this->lower_bandwidth() == m.lower_bandwidth() && this->upper_bandwidth() == m.upper_bandwidth()) {
            // Si tenemos dos matrices banda con los mismos anchos de banda, simplemente sumamos la matriz miembro a miembro.
            int bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;

            for (int i = 0; i < this->rows(); ++i) {
                for (int j = 0; j < bound; ++j) {
                    this->matrix[i][j] += m.matrix[i][j];
                }
            }
        } else {
            // Si no, nos fijamos cuales son los nuevos anchos de banda
            int new_lband = std::max(this->lower_bandwidth(), m.lower_bandwidth());
            int new_uband = std::max(this->upper_bandwidth(), m.lower_bandwidth());
            int new_bound = new_lband + new_uband + 1;

            // Creamos una nueva matriz que guarda directamente la suma
            double **output = new double*[this->rows()];

            for (int i = 0; i < this->rows(); ++i) {
                output[i] = new double[new_bound];

                for (int j = 0; j < new_bound; ++j) {
                    output[i][j] = this->matrix[i][j + i - this->lower_bandwidth()] + m.matrix[i][j + i - m.lower_bandwidth()];
                }
            }

            // Nos aseguramos que los cambios sean efectivos
            this->lband = new_lband;
            this->uband = new_uband;

            // Borramos la matriz vieja
            for (int i = 0; i < this->rows(); ++i) {
                delete[] this->matrix[i];
            }

            delete[] this->matrix;

            // Terminamos de fijar los cambios
            this->matrix = output;
        }
    } else {
        // No podemos sumar
        throw new std::out_of_range("Different dimensions for matrix sum");
    }

    return *this;
}

Matrix &Matrix::operator*=(const double &c) {
    int bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;

    for (int i = 0; i < this->rows(); ++i) {
        for (int j = 0; j < bound; ++j) {
            this->matrix[i][j] *= c;
        }
    }

    return *this;
}

Matrix::~Matrix() {
    for (int i = 0; i < this->rows(); ++i) {
        delete[] this->matrix[i];
    }

    delete[] this->matrix;
}

std::ostream &operator<<(std::ostream &os, const Matrix &m) {
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.columns(); ++j) {
            os << (double)m(i, j) << " ";
        }

        os << std::endl;
    }

    os << std::endl;

    return os;
}

Matrix operator+(const Matrix &m, const Matrix &n) {
    Matrix output(m);
    output += n;
    return output;
}

Matrix operator*(const Matrix &m, const double &c) {
    Matrix output(m);
    output *= c;
    return output;
}

Matrix operator*(const Matrix &m, const Matrix &n) {
    if (m.columns() == n.rows()) {
        int lower = std::max(m.lower_bandwidth(), n.lower_bandwidth());
        int upper = std::max(m.upper_bandwidth(), n.upper_bandwidth());

        Matrix output(m.rows(), n.columns(), lower, upper);

        int diagonal = std::min(output.columns(), output.rows());

        for (int d = 0; d < diagonal; ++d) {
            for (int j = std::max(d - lower, 0); j < std::min(upper, output.columns()); ++j) {
                for (int k = 0; k < output.columns();  ++k) {
                    output(d, j) += m(d, k) * n(k, j);
                }
            }
        }

        return output;
    } else {
        throw new std::out_of_range("Matrix product between incompatible matrices.");
    }
}
