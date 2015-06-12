from os import listdir
from os.path import isfile, join
from os import getcwd
import numpy as np
import math
onlyfiles = [ f for f in listdir(getcwd()) if isfile(join(getcwd(),f)) ]

output = open("output.txt", "w")

for f in onlyfiles:
    if f[-3:len(f)] != "out":
        continue

    fStream = open(f)
    
    values = dict()

    for l in fStream:
        splt = l.split('\t\t\t') # partida en key value
        k = splt[0].strip()

        # partimos ahora los valores
        valores = splt[1].strip().split(' ')
        valores = [int(v) for v in valores]

        mean = np.mean(valores)
        degFred = 1
        if len(valores) <= 1:
            degFred = 0

        std = np.std(valores, ddof =degFred) # n-1 es la muestral
        
        percentile = 1.812 # 95th para t-student (n-1) grados de libertad n = 10
        
        error = percentile * (std/math.sqrt(len(valores)))
        
        values[k] = (mean, error) 
    
    fStream.close()
    output.write(f + "\n")


    #for k,v in values.iteritems():
    #    output.write(k + " " + str(v[0]) + " " + str(v[1]) + "\n")

    output.write("kNN Hit CV ")
    for v in values["kNN Hit CV"]:
        output.write(str((v*100)/4200.0) + " ")
    output.write("\n****************************************************\n")


output.close()

