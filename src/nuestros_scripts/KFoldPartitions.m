function KFoldPartitions(k,observations,outputName)

    c = cvpartition(observations, 'KFold',k)
    fileID = fopen(outputName,'w');

    %iteramos cada fold en i
    for i=1:k
        %para cada fold imprimimos en una linea la particion.

        t = training(c,i)'

        % imprimimos la primera aparte asi podemos ordenar bien las comas
        fprintf(fileID, '%i', t(1));

        for j=2:observations
            fprintf(fileID, ',');
            fprintf(fileID, '%i', t(j));
        end

        fprintf(fileID, '\n');
    end

    fclose(fileID);

end