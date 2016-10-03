function [aliSeq] = createAlignmentMappingTwoProteins(uniprotID1, uniprotID2)

r1 = getUniprotInformation('uniprot',uniprotID1);
r2 = getUniprotInformation('uniprot',uniprotID2);

aliSeq = multialign({r1.sequence, r2.sequence});
    
end