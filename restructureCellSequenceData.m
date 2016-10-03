%Take in an array of sequences, protein names, and rna Names that are cell
%arrays, and return a single matlab struct with fields for Header and
%Sequence
function [tableSeqs] = restructureCellSequenceData(seqArray,protNames,rnaNames)

% Find where there are no sequences available
zEmpty = find(cellfun(@isempty,seqArray));

% Take out those entries from our data
protNames(zEmpty) = [];
rnaNames(zEmpty) = [];
seqArray(zEmpty) = [];

N = length(protNames);
numSeqs = 0;

% Find out how many RNA sequences there are per protein id, and count them
% accordingly
numRNAPerProt = ones(N,1);

for i=1:N
    if iscell(seqArray{i})
        numRNAPerProt(i) = length(seqArray{i});
    end
    numSeqs = numSeqs+numRNAPerProt(i);
end

tableSeqs(numSeqs).Header = '';
tableSeqs(numSeqs).Sequence = '';

idx = 1;
for i=1:N
    
    for j=1:numRNAPerProt(i)
        tableSeqs(idx).Header = char(strcat(protNames{i},'_',rnaNames{i}{j}));
        if iscell(seqArray{i})
            tableSeqs(idx).Sequence = char(seqArray{i}{j});
        else
            tableSeqs(idx).Sequence = char(seqArray{i});
        end
        idx=idx+1;
    end
end

% zEmpty = find(isempty(tableSeqs(:).Sequence));
% tableSeqs(zEmpty).Header = [];
% tableSeqs(zEmpty).Sequence = [];

end