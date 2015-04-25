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

    /*
     * Leemos el training set en una matriz
     */
    std::fstream train("train.csv", std::ios_base::in);

    std::vector<std::vector<char>> matrix;
    while (!train.eof() && train.good()) {
        for (int i = 0; i < DIM*DIM + 1; ++i) {

        }
    }

    if (!train.good()) {
        // TODO: die
    }

    std::fstream test("test.csv", std::ios_base::in);

    std::string path;
    int neighbours, alpha, tests;
    std::fstream input(argv[1], std::ios_base::in);

    input >> path >> neighbours >> alpha >> tests;

    std::bitset<42000> *masks = new std::bitset<42000>[tests]();

    for (int i = 0; i < tests; ++i) {
        input >> masks[i];


        //std::fstream outHandle(output, std::ios_base::out | std::ios_base::ate | std::ios_base::app);

        //if (outHandle.good()) {
        //    outHandle << ;
        //}

    }

    return 0;
}