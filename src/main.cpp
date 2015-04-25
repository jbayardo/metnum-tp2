//
// Created by Julian Bayardo on 4/25/15.
//

#include "Problem.h"
#include <string>
#include <bitset>
#include <fstream>

#define DIM 28

int main(int argc, char *argv[]) {
    SolutionMethod method = SolutionMethod::KNN;

    if (*argv[3] == '1') {
        method = SolutionMethod::PCA_KNN;
    }

    std::string path;
    int neighbours, alpha, tests;
    std::fstream input(argv[1], std::ios_base::in);

    input >> path >> neighbours >> alpha >> tests;

    std::bitset<42000> *masks = new std::bitset<42000>[tests]();

    for (int i = 0; i < tests; ++i) {
        input >> masks[i];
    }

    std::fstream train(path + "train.csv", std::ios_base::in);

    while (!train.eof() && train.good()) {
        for (int i = 0; i < DIM*DIM + 1; ++i) {

        }
    }

    if (!train.good()) {
        // TODO: die
    }

    std::fstream test(path + "test.csv", std::ios_base::in);
    return 0;
}