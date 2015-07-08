import numpy as np
import os
import re
import pandas as pd

# Este script larga output.csv, con el formato que buscamos para que lo tome
# el software para graficar, Tableau.

knnDataset = []

skip = ['Load Testing Dataset Timer', 'Output Dataset Timer', 'Load Training Dataset Timer']
cvPartitions = 10
trainingElements = 4200

def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]

for fname in os.listdir('../resultados'):
    fname = '../resultados/'+fname

    if not os.path.isfile(fname):
        continue

    with open(fname, 'r') as handle:
        data = handle.read()


    matches = re.match(r"(.*)?\/?statistics(knn|pca)_test_([0-9]+)_*([0-9]+)*_*K?([0-9]+)*\.out", fname)
    method = matches.group(2)

    k = int(matches.group(3))

    if method == 'pca':
        alpha = matches.group(4).strip()
    else:
        alpha = None

    K = int(matches.group(5))

    print method, k, alpha, K

    output = {}

    for line in data.split('\n'):
        if len(line) == 0:
            continue

        keyValue = line.split('\t\t\t')

        if keyValue[0] == 'kNN Timer':
            # En este caso tenemos 10 entradas, una por cada
            # entrada en el dataset de training, secuenciales.
            values = reversed(chunks(keyValue[1].strip().split(' '), trainingElements))
        elif (alpha is not None) and (keyValue[0] == 'Power Iteration Iteration Counter' or keyValue[0] == 'Deflation Timer' or keyValue[0] == 'Power Iteration Timer'):
            values = reversed(chunks(keyValue[1].strip().split(' '), int(alpha)))
        elif keyValue[0] not in skip:
            # En este caso tenemos 1 por corrida de cross validation
            values = reversed(keyValue[1].strip().split(' '))
        else:
            continue

        for (index, value) in zip(xrange(cvPartitions), values):
            try:
                output[index]
            except:
                output[index] = {
                    'name': fname,
                    'run': index,
                    'method': method,
                    'k': k,
                    'K': K,
                    'alpha': alpha
                }

            output[index][keyValue[0]] = value

    knnDataset = knnDataset + output.values()

dataset = pd.DataFrame(knnDataset)
dataset.to_csv('tableau.csv', index=False)
