% Takes in the name of a directory (usually a UniProt identifier as found
% in alignments_ecs folder), and returns a Matlab structure with the
% following fields:
% results.uniprot:  The uniprot identifier (should be a cell array)
% results.eval:     an array containing the relevant e-value 
% results.pairs:    a cell array where each entry is the list of EC pairs at
%                   a given e-value
% results.residues: a cell array of the two 1-letter residue codes for
%                   a given e value
% results.ecvals:   a cell array of the EC score returned for each pair,
%                   with each entry corresponding to a different e value

function [results] = readECsIntoMatlab(dirname)

% Navigate to given directory and get list of EC text files
cd(dirname);
listfiles = dir('*_ECs.txt');
listfiles = {listfiles.name};
numFiles = length(listfiles);

if numFiles == 0
    disp(strcat('Error:  No EC files found for ',dirname));
end
results.uniprot = cell(numFiles,1);
results.eval = zeros(numFiles,1);
results.pairs = cell(numFiles,1);
results.residues = cell(numFiles,1);
results.ecvals = cell(numFiles,1);


% Get options that are included in file name (ie gap threshold, etc)
for i=1:numFiles
    tokenNames = regexp(listfiles{i},'(?<uniprot>[^_]+).*_e(?<eval>[0-9]).*','names');
    
    % Store pairs, ecs, etc. in structure
    results.uniprot{i} = tokenNames.uniprot;
    results.eval(i) = str2double(tokenNames.eval);
    [results.pairs{i},results.residues{i},results.ecvals{i}] = readECFile(listfiles{i});
end

end

% Takes in a file containing ECs, sorts them in descending order, and
% stores as a matlab structure
function [pairs,residues,ecvals] = readECFile(inputfile)

indata = readtable(inputfile,'Format','%f%c%f%c%f%f');
indata = sortrows(indata,-6);

pairs = [table2array(indata(:,1)) table2array(indata(:,3))];
residues = [table2array(indata(:,2)) table2array(indata(:,4))];
ecvals = table2array(indata(:,6));

end