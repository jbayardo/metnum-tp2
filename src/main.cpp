//
// Created by Julian Bayardo on 4/25/15.
//

#include <vector>
#include <functional>
#include <bitset>
#include <fstream>
#include <queue>
#include <cmath>
#include "Problem.h"
#include "Matrix.h"

#define DIM 28
#define TRAIN_SIZE 42000
#define TEST_SIZE 28000

typedef unsigned char Label;

void loadTrainingSet(std::string path, Matrix &trainingSet, Label *trainingLabels) {
    std::fstream train(path + "train.csv", std::ios_base::in);
    std::cout << "Levantando training set" << std::endl;
    // Obtenemos el header, así después el algorítmo sólo levanta los datos.
    std::string line;
    std::getline(train, line);

    // Se refiere a la linea del csv que estamos leyendo. O sea, a la imagen que estamos levantando.
    int l = 0;

    // Levantamos una linea del csv y la separamos por coma
    while (getline(train,line) && !train.eof() && !train.bad()) {
        // Buscamos cuál es el Label
        std::string::size_type prev = line.find_first_of(',');

        if (prev == std::string::npos) {
            std::cerr << "Label no encontrado en la linea " << l << " del archivo de training." << std::endl;
            std::cerr << line << std::endl;
            train.close();
            throw std::string("Puto");
        }

        Label lbl = line.substr(0, prev).c_str()[0] - 48;

        if (lbl > 9) {
            std::cerr << "Label erroneo en la linea " << l << " en el archivo de training." << std::endl;
            train.close();
            throw std::string("Culo"); // TODO:
        }

        trainingLabels[l] = lbl;

        // Este contador es el número de pixel que estamos procesando en la imagen
        int i = 0;

        while (prev != std::string::npos) {
            // Encontramos la próxima coma
            std::string::size_type cur = line.find_first_of(',', prev + 1);

            if (cur != std::string::npos) {
                std::string current = line.substr(prev + 1, cur);
                // Levantamos el string como un char
                Label out = (Label) std::stoi(current);
                trainingSet(l, i) = out;
            }

            prev = cur;
            ++i;
        }

        ++l;
    }

    if (train.bad()) {
        std::cerr << "Error de lectura en el archivo de training." << std::endl;
        train.close();
        throw std::string("Pito"); // TODO:
    }

    train.close();
}

void loadTestingSet(std::string path, Matrix &testingSet, Label *testingLabels) {
    std::fstream test(path + "test.csv", std::ios_base::in);
    std::cout << "Levantando testing set" << std::endl;
    // Obtenemos el header, así después el algorítmo sólo levanta los datos.
    std::string line;
    std::getline(test, line);

    // Se refiere a la linea del csv que estamos leyendo. O sea, a la imagen que estamos levantando.
    int l = 0;

    // Levantamos una linea del csv y la separamos por coma
    while (getline(test,line) && !test.eof() && !test.bad()) {
        // Buscamos cuál es el Label
        std::string::size_type prev = line.find_first_of(',');
        testingSet(l, 0) = (Label) std::stoi(line.substr(0, prev));

        // Este contador es el número de pixel que estamos procesando en la imagen
        int i = 1;

        while (prev != std::string::npos) {
            // Encontramos la próxima coma
            std::string::size_type cur = line.find_first_of(',', prev + 1);

            if (cur != std::string::npos) {
                // Levantamos el string como un char
                testingSet(l, i) = (Label) std::stoi(line.substr(prev + 1, cur));
            }

            prev = cur;
            ++i;
        }

        ++l;
    }

    if (test.bad()) {
        std::cerr << "Error de lectura en el archivo de testing." << std::endl;
        test.close();
        throw std::string("Pito"); // TODO:
    }

    test.close();
}

int main(int argc, char *argv[]) {
    SolutionMethod method = SolutionMethod::KNN;

    if (argc < 3) {
        std::cerr << "Faltan argumentos." << std::endl;
        return 0;
    }

    if (*argv[3] == '1') {
        method = SolutionMethod::PCA_KNN;
    }

    // Levantamos los datos de entrada al programa
    std::string path;
    int neighbours, alpha, tests;
    std::fstream input(argv[1], std::ios_base::in);

    input >> path >> neighbours >> alpha >> tests;

    // Levantamos las mascaras para definir el training set y el test set
    std::cout << "Levantando mascaras para el dataset" << std::endl;
    std::bitset<TRAIN_SIZE> *masks = new std::bitset<TRAIN_SIZE>[tests]();

    for (int i = 0; i < tests; ++i) {
        input >> masks[i];
    }

    input.close();

    Matrix trainingSet = Matrix(TRAIN_SIZE, DIM*DIM);
    Label *trainingLabels = new Label[TRAIN_SIZE];
    loadTrainingSet(path, trainingSet, trainingLabels);

    Matrix testingSet = Matrix(TEST_SIZE, DIM*DIM);
    Label *testingLabels = new Label[TEST_SIZE];
    loadTestingSet(path, testingSet, trainingLabels);


    return 0;
}

template <typename T>
using min_queue = std::priority_queue<T, std::vector<T>, std::greater<T>>;

// Matriz de datos, vector, número de línea.
typedef std::function<double(const Matrix &, const std::vector<double> &, int)> DistanceF;

const DistanceF L2 = DistanceF([](const Matrix &A, const std::vector<double> &v, int i) -> double {
    if (v.size() < A.columns()) {
        throw new std::out_of_range("Invalid size for vector");
    }

    double output = 0.0;

    for (int j = 0; j < A.columns(); ++j) {
        output += (static_cast<double>(A(i, j)) - v[i]) * (static_cast<double>(A(i, j)) - v[i]);
    }

    return std::sqrt(output);
});

/*
 * Devuelve un número del 0 al 9 que representa el dígito reconocido
 */
Label kNN(int k, const Matrix &trainingSet, Label *trainingLabels, const std::vector<double> &vector, const DistanceF &f) {
    min_queue<std::pair<double, Label>> distances;

    for (int i = 0; i < trainingSet.rows(); ++i) {
        distances.push(std::pair<double, Label>(f(trainingSet, vector, i), trainingLabels[i]));
    }

    int i = 0;
    int labels[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    while (!distances.empty() && i < k) {
        labels[distances.top().second]++;
        distances.pop();
        ++i;
    }

    Label maximum = 0;

    for (Label j = 0; j < 10; ++j) {
        if (labels[j] > maximum) {
            maximum = j;
        }
    }

    return maximum;
}

typedef std::function<double(const std::vector<double> &v)> Norm;

const Norm N2 = Norm([](const std::vector<double> &v) -> double {
    double output = 0.0;

    for (int i = 0; i < v.size(); ++i) {
        output += i*i;
    }

    return std::sqrt(output);
});

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

/*
 * Obtiene el autovalor dominante en módulo de una matriz.
 *
 * @param A matriz para buscar autoespacios
 * @param x vector inicial del algoritmo
 */
std::pair<double, std::vector<double>> powerIteration(const Matrix &A, const std::vector<double> &x, const Norm &norm, unsigned int iterations = 100) {
    std::vector<double> C(A.columns());

    // Copiamos el vector hasta la coordenada que vamos a usar
    for (int i = 0; i < A.columns(); ++i) {
        C[i] = x[i];
    }

    double eigenValue = 0.0;

    // reutilizamos la norma calculada en el paso anterior.
    double lastNorm = norm(C);

    for (int k = 0; k < iterations; ++k) {
        // Elevamos a potencia
        C = A * C;

        // Normalizamos el vector
        double length = norm(C);        
        for (int i = 0; i < A.columns(); ++i) {
            C[i] /= length;
        }

        eigenValue = length/lastNorm;
        lastNorm = length;
    }

    return std::pair<double, std::vector<double>>(eigenValue, C);
}