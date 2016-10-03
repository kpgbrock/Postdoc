function [newGapTrx,newMrnaIDs,newSeqErrors,chosenTrx] = discriminateAmongTranscripts(cellGapTrx,mrnaIDs,seqErrors)

[N,~] = size(cellGapTrx);
newGapTrx = cell(N,1);
newMrnaIDs = cell(N,1);
newSeqErrors = cell(N,1);
chosenTrx = zeros(N,1);

for i=1:N
    
    % Take out entries where we were not able to return an entry for the
    % EMBL ID
    if ~isempty(cellGapTrx{i})
        
        zEmptySeq = find(cellfun(@isempty,cellGapTrx{i}));
        cellGapTrx{i}(zEmptySeq) = [];
        mrnaIDs{i}(zEmptySeq) = [];
        seqErrors{i}(zEmptySeq) = [];
    end
    
    % If our protein doesn't have any identified transcripts
    if length(seqErrors{i}) <= 1
        newGapTrx{i} = cellGapTrx{i};
        newMrnaIDs{i} = mrnaIDs{i};
        newSeqErrors{i} = seqErrors{i};
        if length(seqErrors{i}) == 1
            chosenTrx(i) = 1;
        end
    else
        [minVal,minI] = min(seqErrors{i});
        newGapTrx{i} = cellGapTrx{i}(minI);
        newMrnaIDs{i} =mrnaIDs{i}(minI);
        newSeqErrors{i} = minVal;
        chosenTrx(i) = minI;
        
    end
    
end

end