function WriteInputFile(outputName, inputPath, k, alpha, partitionCount, cvpart, observations)
 
    fileID = fopen(outputName,'w');
    fprintf(fileID, '%s ', inputPath);
    fprintf(fileID, '%i ', k);
    fprintf(fileID, '%i ', alpha);
    fprintf(fileID, '%i\n', partitionCount);
    
    PrintKFoldPartitions(partitionCount, observations, fileID, cvpart)
    
    fclose(fileID);
end