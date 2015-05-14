//
// Created by Julian Bayardo on 5/8/15.
//

#include <chrono>
#include "LinearAlgebra.h"
#include "Counter.h"

std::vector<double> operator*(const Matrix &m, const std::vector<double> &n) {
    if (n.size() < m.columns()) {
        std::stringstream fmt;
        fmt << "Tamaño de matriz M es " << m.columns() << ", mientras que vector n es " << n.size();
        throw new std::out_of_range(fmt.str());
    }

    std::vector<double> output(m.rows(), 0.0);

    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.columns(); ++j) {
            output[i] += m(i, j) * n[j];
        }
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
EigenPair powerIteration(const Matrix &A, std::vector<double> eigenVector, const Norm &norm, unsigned int iterations) {
    if (A.columns() != A.rows()) {
        throw new std::invalid_argument("La matriz no es cuadrada en el método de la potencia");
    }

    if (iterations <= 0) {
        throw new std::invalid_argument("La cantidad de iteraciones para el método de la potencia debe ser mayor a 0");
    }

    Timer timer("Power Iteration Timer");
    Counter iteration("Power Iteration Iteration Counter");

    double length = norm(eigenVector);

    for (int j = 0; j < eigenVector.size(); ++j) {
        eigenVector[j] /= length;
    }

    // Verificamos convergencia
    while (iteration < iterations) {
        // Elevamos a potencia
        eigenVector = A * eigenVector;

        // Normalizamos el vector
        length = norm(eigenVector);

        for (int j = 0; j < A.rows(); ++j) {
            eigenVector[j] /= length;
        }

        // Actualizamos el contador
        ++iteration;
    }

    double eigenValue = 0.0;
    /*std::vector<double> outfuck(A*eigenVector);

    for (int i = 0; i < outfuck.size(); ++i) {
        eigenValue += outfuck[i] * outfuck[i];
    }*/

    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.rows(); ++j) {
            eigenValue += eigenVector[i] * eigenVector[j] * A(i, j);
        }
    }

    eigenValue /= std::pow(norm(eigenVector), 2);

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
    if (A.columns() != A.rows()) {
        throw new std::runtime_error("La matriz no es cuadrada en el método de deflación");
    }

    Timer timer("Deflation Timer");

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
std::list<EigenPair> decompose(Matrix deflated, int k, const Norm &norm, unsigned int iterations) {
    if (k >= deflated.columns()) {
        std::stringstream fmt;
        fmt << "Cantidad de autovalores esperado es demasiado grande, " << k << " en una matriz de " << deflated.columns();
        throw new std::out_of_range(fmt.str());
    }

    Timer timer("Decompose Timer");
    std::list<EigenPair> output;
    // Vector inicial para esta iteracion
    std::vector<double> x0((unsigned long) deflated.columns(), 0.0);

    for (int i = 0; i < k; ++i) {
        for (int i = 0; i < deflated.columns(); ++i) {
            x0[i] = random() % 1337 + 1;
        }

        // Obtenemos el i-esimo eigenpair dominante
        EigenPair dominant = powerIteration(deflated, x0, norm, iterations);

        if (dominant.first != dominant.first) {
            std::cerr << "Error sacando el autovalor " << i << ". Vector: " << std::endl;

            for (int i = 0; i < deflated.columns(); ++i) {
                std::cerr << x0[i] << " ";
            }

            std::cerr << std::endl;
            --i;
        } else {
            // Hacemos deflacion, para el proximo paso
            deflation(deflated, dominant);

            // Lo guardamos al final de la lista
            output.push_back(dominant);
        }
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

std::vector<double> operator*(const double &m, const std::vector<double> &n) {
    std::vector<double> output(n);
    std::transform(output.begin(), output.end(), output.begin(), [m](const double &x) -> double { return x * m; });
    return output;
}