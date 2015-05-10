//
// Created by Julian Bayardo on 5/8/15.
//

#include "LinearAlgebra.h"

std::vector<double> operator*(const Matrix &m, const std::vector<double> &n) {
    if (n.size() < m.columns()) {
        throw new std::out_of_range("Invalid vector size for matrix product");
    }

    std::vector<double> output = std::vector<double>();

    for (int i = 0; i < m.rows(); ++i) {
        double tmp = 0.0;

        for (int j = 0; j < m.columns(); ++j) {
            tmp += m(i, j) * n[j];
        }

        output[i] = tmp;
    }

    return output;
}

/**
 * Obtiene el autovalor dominante en modulo de una matriz.
 *
 * @param A matriz para buscar autoespacios
 * @param x vector inicial del algoritmo
 * @param norm norma vectorial para actualizar los valores
 * @param condition condicion para verificar la convergencia del metodo
 */
EigenPair powerIteration(const Matrix &A, std::vector<double> eigenVector, const Norm &norm, const ConditionF &condition) {
    // Nos aseguramos que el vector tenga solo las coordenadas que vamos a usar
    eigenVector.resize((unsigned long) A.columns());

    double eigenValue = 0.0;
    int iteration = 0;
    double lastNorm = norm(eigenVector);

    // Verificamos convergencia
    while (!condition(A, eigenVector, eigenValue, iteration)) {
        // Elevamos a potencia
        eigenVector = A * eigenVector;

        // Normalizamos el vector
        double length = norm(eigenVector);

        std::transform(eigenVector.begin(), eigenVector.end(), eigenVector.begin(), [length](const double &x) -> double {
            return x / length;
        });

        // Actualizamos los datos
        eigenValue = length/lastNorm;
        lastNorm = length;

        // Actualizamos el contador
        ++iteration;
    }

    return std::pair<double, std::vector<double>>(eigenValue, eigenVector);
}

/**
 * Corre deflacion sobre la matriz, devuelve una nueva matriz, que retiene los autovectores de la matriz A, excepto por
 * el autovector de autovalor dominante.
 *
 * @param A matriz inicial para la deflacion
 * @param eigen par de autovalor y autovector
 * @return nueva matriz con las caracteristicas anunciadas
 */
void deflation(Matrix &A, const EigenPair &eigen) {
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.columns(); ++j) {
            A(i, j) -= eigen.first * eigen.second[i] * eigen.second[j];
        }
    }
}

/**
 * Obtiene los primeros k autovectores y autovalores dominantes de la matriz.
 *
 * @param A matriz a investigar
 * @param k cantidad de autovectores/autovalores a obtener
 * @ret lista ordenada por dominancia decreciente del autopar
 */
std::list<EigenPair> decompose(Matrix deflated, int k, const Norm &norm, const ConditionF &condition) {
    if (k >= deflated.columns()) {
        throw new std::out_of_range("Too many eigenpairs to obtain");
    }

    std::list<EigenPair> output;

    // Vector inicial para esta iteracion
    std::vector<double> x0((unsigned long) deflated.columns(), 1.0);

    for (int i = 0; i < k; ++i) {
        // Obtenemos el i-esimo eigenpair dominante
        EigenPair dominant = powerIteration(deflated, x0, norm, condition);

        // Hacemos deflacion, para el proximo paso
        deflation(deflated, dominant);

        // Lo guardamos al final de la lista
        output.push_back(dominant);
    }

    return output;
}

/**
 * Realiza un cambio de base a las filas de src, guardandolas en dst y utlizando los autovectores indicados.
 *
 * @param src matriz con imagenes
 * @param dst matriz destino
 * @param l lista de EigenPair
 */

void dimensionReduction(const Matrix& src, Matrix& dst, const std::list<EigenPair>& l)
{
    // contador de EigenPair
    int c = 0;
    // ASUMO QUE LO RECORRE EN ORDEN!
    for (const EigenPair& ep : l)
    {
        const std::vector<double>& eigenVector = ep.second;

        // 
        for (int i = 0; i < src.rows(); i++)
        {   
            double sum = 0.0;
            for (int j = 0; j < src.columns(); j++)
                sum += eigenVector[j] * src(i,j);

            // los guardamos en fila, asi reutilizamos otros metodos.
            dst(i,c) = sum;
        }

        c++;
    }
}