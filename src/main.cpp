//
// Created by Julian Bayardo on 4/25/15.
//

#include <vector>
#include <bitset>
#include <fstream>
#include <sstream>
#include "LinearAlgebra.h"
#include "Counter.h"

#define DIM 28
#define TRAIN_SIZE 42000
#define TEST_SIZE 28000

typedef double Label;

typedef enum {
    KNN,
    PCA_KNN
} SolutionMethod;

/*
 * Devuelve un número del 0 al 9 que representa el dígito reconocido
 */
Label kNN(int k, const Matrix &trainingSet, const std::vector<Label> &trainingLabels, Matrix &evSet, int i1, const DistanceF &f) {
    Counter timer("kNNTimer");
    auto start = std::chrono::steady_clock::now();

    min_queue<std::pair<double, Label>> distances;

    for (int i = 0; i < trainingSet.rows(); ++i) {
        distances.push(std::pair<double, Label>(f(trainingSet, i, evSet, i1), trainingLabels[i]));
    }

    // TODO: considerar posibilidades de clasificacion.
    int i = 0;
    int labels[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    while (!distances.empty() && i < k) {
        labels[(int)distances.top().second]++;
        distances.pop();
        ++i;
    }

    int maximum = 0;

    for (int j = 0; j < 10; ++j) {
        if (labels[j] > maximum) {
            maximum = j;
        }
    }

    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    timer.set(elapsed.count());

    return maximum;
}

void loadTrainingSet(std::string path, Matrix &trainingSet, std::vector<Label> &trainingLabels) {
    Counter timer("loadTrainingSetTimer");
    auto start = std::chrono::steady_clock::now();

    std::fstream train(path + "train.csv", std::ios_base::in);
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
            std::stringstream fmt;
            fmt << "Label no encontrado en la linea " << l << " en el archivo de training";
            train.close();
            throw new std::invalid_argument(fmt.str());
        }

        Label lbl = line.substr(0, prev).c_str()[0] - 48;

        if (lbl > 9) {
            std::stringstream fmt;
            fmt << "Label erroneo en la linea " << l << " en el archivo de training";
            train.close();
            throw new std::out_of_range(fmt.str());
        }

        trainingLabels[l] = lbl;

        // Este contador es el número de pixel que estamos procesando en la imagen
        int i = 0;

        while (prev != std::string::npos) {
            // Encontramos la próxima coma
            std::string::size_type cur = line.find_first_of(',', prev + 1);

            if (cur != std::string::npos) {
                trainingSet(l, i) = std::stod(line.substr(prev + 1, cur));
            }

            prev = cur;
            ++i;
        }

        ++l;
    }

    if (train.bad()) {
        train.close();
        throw new std::runtime_error("Error de lectura en el archivo de training");
    }

    train.close();

    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    timer.set(elapsed.count());
}

void loadTestingSet(std::string path, Matrix &testingSet) {
    Counter timer("loadTestingSetTimer");
    auto start = std::chrono::steady_clock::now();

    std::fstream test(path + "test.csv", std::ios_base::in);
    // Obtenemos el header, así después el algorítmo sólo levanta los datos.
    std::string line;
    std::getline(test, line);

    // Se refiere a la linea del csv que estamos leyendo. O sea, a la imagen que estamos levantando.
    int l = 0;

    // Levantamos una linea del csv y la separamos por coma
    while (getline(test,line) && !test.eof() && !test.bad()) {
        // Buscamos cuál es el primer pixel
        std::string::size_type prev = line.find_first_of(',');
        testingSet(l, 0) = (double)std::stoi(line.substr(0, prev));

        // Este contador es el número de pixel que estamos procesando en la imagen
        int i = 1;

        while (prev != std::string::npos) {
            // Encontramos la próxima coma
            std::string::size_type cur = line.find_first_of(',', prev + 1);

            if (cur != std::string::npos) {
                testingSet(l, i) = std::stod(line.substr(prev + 1, cur));
            }

            prev = cur;
            ++i;
        }

        ++l;
    }

    if (test.bad()) {
        test.close();
        throw new std::runtime_error("Error de lectura en el archivo de testing");
    }

    test.close();

    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    timer.set(elapsed.count());
}

template <std::size_t K>
std::pair<Matrix, std::vector<Label>> filterDataset(const Matrix &A, const std::vector<Label> &labels, const std::bitset<K> &filter) {
    Counter timer("filterDatasetTimer");
    auto start = std::chrono::steady_clock::now();

    if (K < labels.size()) {
        std::stringstream fmt;
        fmt << "Filtro para el dataset tiene tamaño " << K << " cuando el dataset tiene tamaño " << labels.size();
        throw new std::invalid_argument(fmt.str());
    }

    std::vector<Label> output(filter.count(), 0.0);
    int last = 0;

    for (int i = 0; i < A.rows(); ++i) {
        if (filter[i]) {
            output[last] = labels[i];
            last++;
        }

        if (last > filter.count()) {
            break;
        }
    }

    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    timer.set(elapsed.count());

    return std::pair<Matrix, std::vector<Label>>(Matrix(A, filter), output);
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

    std::cerr << "Input file: " << argv[1] << std::endl;
    std::cerr << "Output file: " << argv[2] << std::endl;

    // Levantamos los datos de entrada al programa
    std::string path;
    int neighbours, alpha, tests;
    std::fstream input(argv[1], std::ios_base::in);

    input >> path >> neighbours >> alpha >> tests;

    std::cerr << "Input path: " << path << std::endl;
    std::cerr << "Neighbours: " << neighbours << std::endl;
    std::cerr << "Alpha: " << alpha << std::endl;
    std::cerr << "Tests: " << tests << std::endl;

    // Levantamos las mascaras para definir el training set y el test set
    std::vector<std::bitset<TRAIN_SIZE>> masks((unsigned long) tests);

    for (int i = 0; i < tests; ++i) {
        input >> masks[i];
    }

    input.close();

    Matrix trainingSet(TRAIN_SIZE, DIM*DIM);
    std::vector<Label> trainingLabels(TRAIN_SIZE, 0.0);
    loadTrainingSet(path, trainingSet, trainingLabels);

    Matrix testingSet(TEST_SIZE, DIM*DIM);
    loadTestingSet(path, testingSet);

    std::vector<Label> predictions(TEST_SIZE, 0.0);

    switch (method) {
        case PCA_KNN:
        case KNN:
            std::cerr << "Comenzando kNN para las particiones" << std::endl;

            tests = 0; // TODO:

            for (int k = 0; k < tests; ++k) {
                std::cerr << "Particion " << k << std::endl;
                std::pair<Matrix, std::vector<Label>> fTrain = filterDataset(trainingSet, trainingLabels, masks[k]);
                masks[k].flip();
                std::pair<Matrix, std::vector<Label>> fTest = filterDataset(trainingSet, trainingLabels, masks[k]);
                masks[k].flip();

                Counter hit("kNNHit");
                Counter miss("kNNMiss");

                for (int i = 0; i < fTest.first.rows(); ++i) {
                    Label l = kNN(neighbours, fTrain.first, fTrain.second, fTest.first, i, L2);

                    if (l == fTest.second[i]) {
                        ++hit;
                    } else {
                        ++miss;
                    }
                }
            }

            std::cerr << "Comenzando kNN para el testing set" << std::endl;

            for (int i = 0; i < testingSet.rows(); ++i) {
                Label l = kNN(neighbours, trainingSet, trainingLabels, testingSet, i, L2);
                predictions[i] = l;

                std::cerr << "Predicción " << i << ": " << l << std::endl;

                if (i % 100 == 0) {
                    std::cerr << "Progreso: " << i << std::endl;
                }
            }

            break;
    }

    Counter timer("outputDatasetTimer");
    auto start = std::chrono::steady_clock::now();

    std::fstream output(std::string(argv[2]) + ".csv", std::ios_base::out);
    output << "ImageId,Label" << std::endl;

    for (int i = 0; i < predictions.size(); ++i) {
        output << i << "," << predictions[i] << std::endl;
    }

    output.close();

    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    timer.set(elapsed.count());

    Logger::getInstance().dump("times.txt");

    return 0;
}
