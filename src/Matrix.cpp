//
// Created by Julian Bayardo on 4/25/15.
//

#include "Matrix.h"

const unsigned char zero = 0;

Matrix::Matrix(const Matrix &m) : N(m.rows()), M(m.columns()), uband(m.upper_bandwidth()), lband(m.lower_bandwidth()) {
    int bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;
    this->matrix = new unsigned char*[this->rows()];

    for (int i = 0; i < this->rows(); ++i) {
        this->matrix[i] = new unsigned char[bound];

        for (int j = 0; j < bound; ++j) {
            this->matrix[i][j] = m.matrix[i][j];
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
    this->matrix = new unsigned char*[this->rows()];

    for (int i = 0; i < this->rows(); ++i) {
        this->matrix[i] = new unsigned char[bound];

        for (int j = 0; j < bound; ++j) {
            this->matrix[i][j] = 0.0;
        }
    }
}

int inline Matrix::rows() const {
    return this->N;
}

int inline Matrix::columns() const {
    return this->M;
}

int inline Matrix::upper_bandwidth() const {
    return this->uband;
}

int inline Matrix::lower_bandwidth() const {
    return this->lband;
}

unsigned char &Matrix::operator()(const int &i, const int &j) {
    if (i >= this->rows() || j >= this->columns()) {
        throw new std::out_of_range("Index access out of range");
    }

    if (i <= j + this->lower_bandwidth() && j <= i + this->upper_bandwidth()) {
        return matrix[i][j - i + this->lower_bandwidth()];
    } else {
        throw new std::out_of_range("Out of modifiable range");
    }
}

const unsigned char &Matrix::operator()(const int &i, const int &j) const {
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
        this->matrix = new unsigned char*[m.rows()];

        for (int i = 0; i < this->rows(); ++i) {
            this->matrix[i] = new unsigned char[bound];

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
            unsigned char **output = new unsigned char*[this->rows()];

            for (int i = 0; i < this->rows(); ++i) {
                output[i] = new unsigned char[new_bound];

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

Matrix &Matrix::operator*=(const unsigned char &c) {
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
            os << (int)m(i, j) << " ";
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

Matrix operator*(const Matrix &m, const unsigned char &c) {
    Matrix output(m);
    output *= c;
    return output;
}

Matrix operator*(const Matrix &m, const Matrix &n) {
    if (m.columns() == n.rows()) {
        int max_lower_upper = std::max(m.lower_bandwidth(), n.upper_bandwidth());
        int max_upper_lower = std::max(m.upper_bandwidth(), n.lower_bandwidth());

        Matrix output(m.rows(), n.columns(), max_lower_upper, max_upper_lower);

        int diagonal = std::min(output.columns(), output.rows());

        for (int d = 0; d < diagonal; ++d) {
            for (int j = d - output.lower_bandwidth(); j < std::min(d + output.upper_bandwidth(), output.columns()); ++j) {

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