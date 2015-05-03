//
// Created by Julian Bayardo on 4/25/15.
//

#include "Problem.h"
#include "Matrix.h"
#include <bitset>
#include <fstream>
#include <functional>
#include <queue>

#define DIM 28
#define TRAIN_SIZE 42000
#define TEST_SIZE 28000
// Matriz de datos, vector, número de línea.
typedef std::function<double(const Matrix &, const Matrix &, int i)> Norm;
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

/*
 * Devuelve un número del 0 al 9 que representa el dígito reconocido
 */
Label kNN(int k, const Matrix &trainingSet, Label *trainingLabels, const Matrix &vector, Norm &f) {
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
