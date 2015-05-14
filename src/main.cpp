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
    Timer timer("kNN Timer");
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

    return maximum;
}

void loadTrainingSet(std::string path, Matrix &trainingSet, std::vector<Label> &trainingLabels) {
    Timer timer("Load Training Dataset Timer");

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
            fmt << "Label erroneo en la linea " << l << " en el archivo de training: " << lbl;
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
}

void loadTestingSet(std::string path, Matrix &testingSet) {
    Timer timer("Load Testing Dataset Timer");

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
        testingSet(l, 0) = std::stod(line.substr(0, prev));

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
}

template <std::size_t K>
std::pair<Matrix, std::vector<Label>> filterDataset(const Matrix &A, const std::vector<Label> &labels, const std::bitset<K> &filter) {
    Timer timer("Filter Dataset Timer");

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

    return std::pair<Matrix, std::vector<Label>>(Matrix(A, filter), output);
}

void PCAKNN(std::string path, std::string output, std::string append, int alpha, int neighbours, int tests, std::vector<std::bitset<TRAIN_SIZE>> &masks,
            Matrix &trainingSet, std::vector<Label> &trainingLabels, Matrix &testingSet, std::vector<Label> &predictions) {
    std::cerr << "Comenzando kNN con PCA para las particiones" << std::endl;

    for (int k = 0; k < tests; ++k) {
        std::cerr << "Procesando partición " << k << std::endl;

        std::pair<Matrix, std::vector<Label>> fTrain = filterDataset(trainingSet, trainingLabels, masks[k]);
        std::cerr << "Casos de train: " << masks[k].count() << std::endl;
        std::cerr << "Casos de train concretos: " << fTrain.first.rows() << std::endl;

        masks[k].flip();
        std::pair<Matrix, std::vector<Label>> fTest = filterDataset(trainingSet, trainingLabels, masks[k]);
        std::cerr << "Casos de test: " << masks[k].count() << std::endl;
        std::cerr << "Casos de test concretos: " << fTest.first.rows() << std::endl;

        Timer timer("kNN Partition Timer");

        // Las direcciones en las que hay mayor dispersión de datos son los autovectores de la matriz de covarianza.
        // Estos autovectores forman una base ortonormal.
        // Vamos a obtenerlos (cierta cantidad) y luego realizar cambio de base, para realizar KNN en menos dimensiones.

        // Debemos armar la matriz de covarianza, de las imagenes de training test.

        // primero debemos calcular el promedio de nuestras imagenes.
        // lo definimos como una fila.
        std::cerr << "Calculando promedio de variables para la particion " << k << std::endl;
        Matrix mean(1, DIM*DIM);
        Timer PCAMean("PCA Mean Testing Dataset");

        // TODO: Seria util que la propia matriz tenga una funcion que te de el vector promedio. La sumatoria la va generando cada vez que se cambian las filas.
        for (int j = 0; j < fTrain.first.columns(); j++) {
            for (int i = 0; i < fTrain.first.rows(); i++) {
                mean(0, j) += fTrain.first(i, j);
            }

            mean(0, j) /= fTrain.first.rows();
        }

        PCAMean.stop();

        Matrix covariance(fTrain.first.columns(), fTrain.first.columns());
        std::string covFileName = path + "covariance" + std::to_string(k) + append; // TODO: convertir k a string
        std::fstream inCov(covFileName, std::ios_base::in);

        if (inCov.good()) {
            std::cerr << "Levantando matriz de covarianza para la partición " << k << " del training dataset" << std::endl;
            inCov >> covariance;
            inCov.close();
        } else {
            std::cerr << "Calculando matriz de covarianza la partición " << k << " training dataset" << std::endl;
            Timer PCACov("PCA Covariance Testing Dataset");

            // TODO: hacemos dos veces la diagonal
            for (int j = 0; j < covariance.columns(); j++) {
                if (j % 10 == 0) {
                    std::cerr << "Progreso: " << j << "/" << covariance.columns() << std::endl;
                }

                double m1 = mean(0, j);

                for (int k = 0; k <= j; k++) {
                    double m2 = mean(0, k);

                    for (int i = 0; i < fTrain.first.rows(); i++) {
                        covariance(j, k) += (fTrain.first(i, j) - m1) * (fTrain.first(i, k) - m2);
                    }

                    covariance(j, k) /= (fTrain.first.rows() - 1 );
                    covariance(k, j) = covariance(j, k);
                }
            }

            PCACov.stop();

            std::fstream outCov(covFileName, std::ios_base::out);

            if (outCov.good()) {
                std::cerr << "Guardando matriz de covarianza para la partición " << k << " del training dataset" << std::endl;
                outCov << covariance;
                outCov.close();
            } else {
                std::cerr << "Error guardando la matriz de covarianza, siguiendo igualmente" << std::endl;
            }
        }

        // no nos importa buscar los autovalores, solo los autovectores.

        std::list<EigenPair> eigenPair = decompose(covariance, alpha, N2, 1000);
        // asumimos que aca tenemos los autovectores y autovalores

        std::string eigenFileName = path + "eigenvalues" + std::to_string(k) + append; // TODO: convertir k a string
        std::fstream eigenvalues(eigenFileName, std::ios_base::out);

        if (eigenvalues.good()) {
            eigenvalues.precision(7);
            std::cerr << "Guardando autovalores para los test." << std::endl;

            for (const EigenPair& ep : eigenPair) {
                eigenvalues << ep.first << std::endl;
            }

            eigenvalues.close();
        } else {
            std::cerr << "Error guardando los autovalores, siguiendo igualmente" << std::endl;
        }

        std::cerr << "Haciendo cambio de base para el training dataset" << std::endl;
        Matrix trainChangeBasis(fTrain.first.rows(), alpha);
        // En este paso vamos a realizar un cambio de espacio a todos los vectores
        dimensionReduction(fTrain.first, trainChangeBasis, eigenPair);

        // a cada imagen del testing set debemos restarle mean y dividirlos por sqrt(fTrain.first.rows()-1)
        // segun diapositivas de la clase.

        for (int i = 0; i < fTest.first.rows(); i++) {
            for (int j = 0; j < fTest.first.columns(); j++) {
                fTest.first(i,j) -= mean(0, j);
                fTest.first(i,j) /= sqrt(fTrain.first.rows()-1);
            }
        }

        std::cerr << "Haciendo cambio de base para el testing dataset" << std::endl;
        Matrix testChangeBasis(fTest.first.rows(), alpha);
        dimensionReduction(fTest.first, testChangeBasis, eigenPair);

        // ya tenemos los vectores en sus respectivos cambios de bases
        Counter hit("kNN Hit");
        Counter miss("kNN Miss");

        for (int i = 0; i < testChangeBasis.rows(); ++i) {
            Label l = kNN(neighbours, trainChangeBasis, fTest.second, testChangeBasis, i, L2);

            if (l == fTest.second[i]) {
                ++hit;
            } else {
                ++miss;
            }
        }
    }

    std::cerr << "Comenzando PCA para el testing dataset" << std::endl;
    // Las direcciones en las que hay mayor dispersión de datos son los autovectores de la matriz de covarianza.
    // Estos autovectores forman una base ortonormal.
    // Vamos a obtenerlos (cierta cantidad) y luego realizar cambio de base, para realizar KNN en menos dimensiones.

    // Debemos armar la matriz de covarianza, de las imagenes de training test.

    // primero debemos calcular el promedio de nuestras imagenes.
    // lo definimos como una fila.
    std::cerr << "Calculando promedio de variables para el training dataset" << std::endl;
    Matrix mean(1, DIM*DIM);
    Timer PCAMean("PCA Mean Testing Dataset");

    // TODO: Seria util que la propia matriz tenga una funcion que te de el vector promedio. La sumatoria la va generando cada vez que se cambian las filas.
    for (int i = 0; i < trainingSet.rows(); i++) {
        for (int j = 0; j < trainingSet.columns(); j++) {
            mean(0,j) += trainingSet(i,j)/TRAIN_SIZE;
        }
    }

    PCAMean.stop();

    std::cerr << "Normalizando training dataset" << std::endl;
    Timer PCASum("PCA Sum Testing Dataset");

    // Debemos generar en la matriz de training lo siguiente en cada fila:
    // x_i debe ser (x_i - mean)_traspuesto / (sqrt(n-1))
    // para nuestro caso ya estan traspuestas.
    for (int i = 0; i < trainingSet.rows(); i++) {
        for (int j = 0; j < trainingSet.columns(); j++) {
            trainingSet(i,j) -= mean(0, j);
            trainingSet(i,j) /= sqrt(TRAIN_SIZE-1);
        }
    }

    PCASum.stop();

    Matrix covariance(trainingSet.columns(), trainingSet.columns());

    std::fstream inCov(path + "covariance.ssv", std::ios_base::in);

    if (inCov.good()) {
        std::cerr << "Levantando matriz de covarianza para el training dataset" << std::endl;
        inCov >> covariance;
        inCov.close();
    } else {
        std::cerr << "Calculando matriz de covarianza para el training dataset" << std::endl;
        Timer PCACov("PCA Covariance Testing Dataset");

        // TODO: hacemos dos veces la diagonal
        for (int i = 0; i < covariance.columns(); i++) {
            if (i % 10 == 0) {
                std::cout << "Progreso: " << i << std::endl;
            }

            for (int j = 0; j <= i; j++) {
                // j es la columna de X_t, que resulta ser la fila j-esima de X
                for (int k = 0; k < trainingSet.rows(); k++) {
                    covariance(i, j) += trainingSet(k, j) * trainingSet(k, i);
                    covariance(j, i) = covariance(i, j);
                }
            }
        }

        PCACov.stop();

        std::fstream outCov(path + "covariance.ssv", std::ios_base::out);

        if (outCov.good()) {
            std::cerr << "Guardando matriz de covarianza para el training dataset" << std::endl;
            outCov << covariance;
            outCov.close();
        } else {
            std::cerr << "Error guardando la matriz de covarianza, siguiendo igualmente" << std::endl;
        }
    }

    // no nos importa buscar los autovalores, solo los autovectores.

    std::list<EigenPair> eigenPair = decompose(covariance, alpha, N2, 100);
    // asumimos que aca tenemos los autovectores y autovalores

    std::fstream eigenvalues(output, std::ios_base::out);

    if (eigenvalues.good()) {
        std::cerr << "Guardando autovalores para los test." << std::endl;

        for (const EigenPair& ep : eigenPair) {
            eigenvalues << ep.first << std::endl;
        }

        eigenvalues.close();
    } else {
        std::cerr << "Error guardando los autovalores, siguiendo igualmente" << std::endl;
    }

    std::cerr << "Haciendo cambio de base para el training dataset" << std::endl;
    Matrix trainChangeBasis(TRAIN_SIZE, alpha);
    // En este paso vamos a realizar un cambio de espacio a todos los vectores
    dimensionReduction(trainingSet, trainChangeBasis, eigenPair);

    // a cada imagen del testing set debemos restarle mean y dividirlos por sqrt(TRAIN_SIZE-1)
    // segun diapositivas de la clase.

    for (int i = 0; i < testingSet.rows(); i++) {
        for (int j = 0; j < testingSet.columns(); j++) {
            testingSet(i,j) -= mean(0, j);
            testingSet(i,j) /= sqrt(TRAIN_SIZE-1);
        }
    }

    std::cerr << "Haciendo cambio de base para el testing dataset" << std::endl;
    Matrix testChangeBasis(TEST_SIZE, alpha);
    dimensionReduction(testingSet, testChangeBasis, eigenPair);

    // ya tenemos los vectores en sus respectivos cambios de bases
    for (int i = 0; i < testChangeBasis.rows(); ++i) {
        Label l = kNN(neighbours, trainChangeBasis, trainingLabels, testChangeBasis, i, L2);
        predictions[i] = l;
    }
}

void NORMALKNN(int alpha, int neighbours, int tests, std::vector<std::bitset<TRAIN_SIZE>> &masks,
               Matrix &trainingSet, std::vector<Label> &trainingLabels, Matrix &testingSet, std::vector<Label> &predictions) {
    std::cerr << "Comenzando kNN para las particiones" << std::endl;

    for (int k = 0; k < tests; ++k) {
        std::cerr << "Procesando partición " << k << std::endl;
        std::pair<Matrix, std::vector<Label>> fTrain = filterDataset(trainingSet, trainingLabels, masks[k]);
        masks[k].flip();
        std::pair<Matrix, std::vector<Label>> fTest = filterDataset(trainingSet, trainingLabels, masks[k]);

        Counter hit("kNN Hit");
        Counter miss("kNN Miss");
        Timer timer("kNN Partition Timer");

        for (int i = 0; i < fTest.first.rows(); ++i) {
            Label l = kNN(neighbours, fTrain.first, fTrain.second, fTest.first, i, L2);

            if (l == fTest.second[i]) {
                ++hit;
            } else {
                ++miss;
            }
        }
    }

    std::cerr << "Comenzando kNN para el testing dataset" << std::endl;

    Timer timer("kNN Testing Timer");

    for (int i = 0; i < testingSet.rows(); ++i) {
        Label l = kNN(neighbours, trainingSet, trainingLabels, testingSet, i, L2);
        predictions[i] = l;

        if (i % 100 == 0) {
            std::cerr << "Progreso de kNN sobre el testing dataset: " << i << "/" << testingSet.rows() << std::endl;
        }
    }
}

std::string basename(std::string const& pathname) {
    return std::string(std::find_if( pathname.rbegin(), pathname.rend(), [](char ch) -> bool { return ch == '/'; }).base(), pathname.end());
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
    std::cerr << "Basename: " << basename(std::string(argv[2])) << std::endl;

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
        for (int j = 0; j < TRAIN_SIZE; ++j) {
            bool in = false;
            input >> in;
            masks[i][j] = in;
        }
    }

    input.close();

    Matrix trainingSet(TRAIN_SIZE, DIM*DIM);
    std::vector<Label> trainingLabels(TRAIN_SIZE, 0.0);
    loadTrainingSet(path, trainingSet, trainingLabels);

    Matrix testingSet(TEST_SIZE, DIM*DIM);
    loadTestingSet(path, testingSet);

    std::vector<Label> predictions(TEST_SIZE, 0.0);

    if (method == PCA_KNN) {
        PCAKNN(path, std::string(argv[2]), basename(std::string(argv[2])), alpha, neighbours, tests, masks, trainingSet, trainingLabels, testingSet, predictions);
    } else {
        NORMALKNN(alpha, neighbours, tests, masks, trainingSet, trainingLabels, testingSet, predictions);
    }

    Timer timer("Output Dataset Timer");

    std::fstream output(path + "clasification.csv", std::ios_base::out);
    output << "ImageId,Label" << std::endl;

    for (int i = 0; i < predictions.size(); ++i) {
        output << i + 1 << "," << predictions[i] << std::endl;
    }

    output.close();
    timer.stop();

    Logger::getInstance().dump("times.txt");

    return 0;
}
