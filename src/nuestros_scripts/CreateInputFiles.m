function CreateInputFiles(neighbours, outputName, inputPath, alpha, partitionCount, observations)
    % En esta funcion generamos cada archivo para cada k dentro de
    % neighbours, usando las mismas particiones.
    % Es decir: dados distintos k, generamos casos de pruebas para las mismas
    % particiones
    
    c = cvpartition(observations, 'KFold', partitionCount);
    for i=neighbours
        WriteInputFile(strcat(outputName, '_', int2str(i), '.in'), inputPath, i, alpha, partitionCount, c, observations);
    end
end
