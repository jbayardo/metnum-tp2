//
// Created by Julian Bayardo on 4/25/15.
//

#include <vector>
#include <bitset>
#include <fstream>
#include "LinearAlgebra.h"

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
                Label out = (double)std::stoi(current);
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
        throw new std::string("Pito"); // TODO:
    }

    train.close();
}

void loadTestingSet(std::string path, Matrix &testingSet) {
    std::fstream test(path + "test.csv", std::ios_base::in);
    std::cout << "Levantando testing set" << std::endl;
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
                // Levantamos el string como un double
                testingSet(l, i) = (double) std::stoi(line.substr(prev + 1, cur));
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

template <std::size_t K>
std::pair<Matrix, std::vector<Label>> filterDataset(const Matrix &A, const std::vector<Label> &labels, const std::bitset<K> &filter) {
    if (K < labels.size()) {
        throw new std::string("Invalid filter for dataset");
    }

    std::vector<Label> output(filter.count());
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
    std::vector<std::bitset<TRAIN_SIZE>> masks((unsigned long) tests);

    for (int i = 0; i < tests; ++i) {
        input >> masks[i];

        // Si el test esta al pedo, nos lo salteamos
        if (masks[i].count() == masks[i].size() || masks[i].count() == 0) {
            --i;
            --tests;
        }
    }

    input.close();

    Matrix trainingSet = Matrix(TRAIN_SIZE, DIM*DIM);
    std::vector<Label> trainingLabels(TRAIN_SIZE);
    loadTrainingSet(path, trainingSet, trainingLabels);

    Matrix testingSet = Matrix(TEST_SIZE, DIM*DIM);
    loadTestingSet(path, testingSet);

    std::vector<Label> predictions;

    switch (method) {
        case PCA_KNN:
            {
                // TODO: CROSS VALIDATION 

                // Las direcciones en las que hay mayor dispersión de datos son los autovectores de la matriz de covarianza.
                // Estos autovectores forman una base ortonormal.
                // Vamos a obtenerlos (cierta cantidad) y luego realizar cambio de base, para realizar KNN en menos dimensiones.

                // Debemos armar la matriz de covarianza, de las imagenes de training test.

                // primero debemos calcular el promedio de nuestras imagenes.
                // lo definimos como una fila.
                Matrix mean(1, DIM*DIM);
                // VERIFICAR!!! LA PRIMERA FILA ES LA 0 ?

                // TODO: Seria util que la propia matriz tenga una funcion que te de el vector promedio. La sumatoria la va generando cada vez que se cambian las filas.
                for (int i = 0; i < trainingSet.rows(); i++)
                    for (int j = 0; j < trainingSet.columns(); j++)
                        mean(0,j) = mean(0,j) + trainingSet(i,j)/TRAIN_SIZE;

                // Debemos generar en la matriz de training lo siguiente en cada fila:
                // x_i debe ser (x_i - mean)_traspuesto / (sqrt(n-1))
                // para nuestro caso ya estan traspuestas.
                for (int i = 0; i < trainingSet.rows(); i++)
                {
                    for (int j = 0; j < trainingSet.columns(); j++)
                    {
                       trainingSet(i,j) -= mean(0, j);
                       trainingSet(i,j) /= sqrt(TRAIN_SIZE-1);
                    }
                }

                Matrix covariance(TRAIN_SIZE, TRAIN_SIZE);
                for (int j = 0; j < TRAIN_SIZE; j++)
                {
                    for (int i = 0; i < TRAIN_SIZE; i++)
                    {  
                        // j es la columna de X_t, que resulta ser la fila j-esima de X
                        for (int k = 0; k < DIM*DIM; k++)
                            covariance(j,i) += trainingSet(k,j) * trainingSet(k,i);
                    }
                }

                // la matriz covariance es simetrica
                // no nos importa buscar los autovalores, solo los autovectores.

                // FIX: No puedo llamarla bien, ni con funciones placeholders.
                //std::list<EigenPair> eigenPair = decompose(covariance, alpha, , );
                std::list<EigenPair> eigenPair; // ASIGNARLES VALORES!
                // asumimos que aca tenemos los autovectores y autovalores


                Matrix trainChangeBasis(TRAIN_SIZE, alpha);
                // En este paso vamos a realizar un cambio de espacio a todos los vectores
                dimensionReduction(trainingSet, trainChangeBasis, eigenPair);

                // a cada imagen del testing set debemos restarle mean y dividirlos por sqrt(TRAIN_SIZE-1)
                // segun diapositivas de la clase.

                for (int i = 0; i < testingSet.rows(); i++)
                {
                    for (int j = 0; j < testingSet.columns(); j++)
                    {
                        testingSet(i,j) -= mean(0, j);
                        testingSet(i,j) /= sqrt(TRAIN_SIZE-1);
                    }
                }

                Matrix testChangeBasis(TEST_SIZE, alpha);
                dimensionReduction(testingSet, testChangeBasis, eigenPair);
                
                // ya tenemos los vectores en sus respectivos cambios de bases
                for (int i = 0; i < testChangeBasis.rows(); ++i) {
                    Label l = kNN(neighbours, trainChangeBasis, trainingLabels, testChangeBasis, i, L2);
                    predictions[i] = l;
                }
                break;
            }
        case KNN:
            for (int k = 0; k < tests; ++k) {
                std::pair<Matrix, std::vector<Label>> fTrain = filterDataset(trainingSet, trainingLabels, masks[k]);
                masks[k].flip();
                std::pair<Matrix, std::vector<Label>> fTest = filterDataset(trainingSet, trainingLabels, masks[k]);
                masks[k].flip();
                unsigned int hit = 0;
                unsigned int miss = 0;

                for (int i = 0; i < fTest.first.rows(); ++i) {
                    Label l = kNN(neighbours, fTrain.first, fTrain.second, fTest.first, i, L2);

                    if (l == fTest.second[i]) {
                        ++hit;
                    } else {
                        ++miss;
                    }
                }

                // TODO: cuentitas.
            }

            for (int i = 0; i < testingSet.rows(); ++i) {
                Label l = kNN(neighbours, trainingSet, trainingLabels, testingSet, i, L2);

                predictions[i] = l;
            }

            break;
    }

    // TODO: output csv


    std::fstream output(std::string(argv[2]) + ".csv", std::ios_base::out);
    output << "ImageId,Label" << std::endl;

    for (int i = 0; i < predictions.size(); ++i) {
        output << i << "," << predictions[i] << std::endl;
    }

    output.close();

    return 0;
}
