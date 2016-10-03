function [seq5,seq3] = get5and3UTR(cellSeq,cellCDS,chosenIndex)

N = length(cellSeq);
seq5 = cell(N,1);
seq3 = cell(N,1);

for i=1:N
    
    %chosenIndex(i)
    if chosenIndex(i) ~= 0
        cds = cellCDS{i}(chosenIndex(i),:);
        
        
        seq5{i} = cellSeq{i}{chosenIndex(i)}(1:cds(1)-1);
        seq3{i} = cellSeq{i}{chosenIndex(i)}(cds(1,2)+1:end);
    end
end
end