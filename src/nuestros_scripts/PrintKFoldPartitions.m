function PrintKFoldPartitions(k,observations,fileID, cvpart)
    
    %iteramos cada fold en i
    for i=1:k
        %para cada fold imprimimos en una linea la particion.

        t = training(cvpart,i)'

        % imprimimos la primera aparte asi podemos ordenar bien las comas
        fprintf(fileID, '%i', t(1));

        for j=2:observations
            fprintf(fileID, ' ');
            fprintf(fileID, '%i', t(j));
        end

        fprintf(fileID, '\n');
    end
end