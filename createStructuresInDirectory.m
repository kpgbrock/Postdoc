function [results] = createStructuresInDirectory(dirName)

currWD = pwd;
addpath(currWD);

cd(dirName);
listFolders = dir();
dirFlags = [listFolders.isdir];
listFolders = {listFolders(dirFlags).name};
listFolders = listFolders(3:end);

numFolders = length(listFolders);

results.protein = cell(numFolders,1);
results.data = cell(numFolders,1);

for i=1:numFolders
    results.protein{i} = listFolders{i};
    results.data{i} = readECsIntoMatlab(listFolders{i});
end

cd(currWD);
end