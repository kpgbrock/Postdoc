function [results] = getmRNAFromProteinMSA(alnFile)

proteinMSA = fastaread(alnFile);
numProts = length(proteinMSA);

results.uniprotIDs = cell(numProts,1);
results.uniprotData = cell(numProts,1);
results.listOfMRNAs = cell(numProts,1);
results.mrnaSeqs = cell(numProts,1);
results.mrnaCDS = cell(numProts,1);
results.mrnaSeqGapped = cell(numProts,1);
results.mrnaSeqErrors = cell(numProts,1);
results.mrnaSeqMismatch = cell(numProts,1);
results.listOfGenDNAs = cell(numProts,1);
results.genSeqs = cell(numProts,1);
results.genCDS = cell(numProts,1);
results.genSeqGapped = cell(numProts,1);
results.genSeqErrors = cell(numProts,1);
results.genSeqMismatch = cell(numProts,1);

% Need to set up first protein, our query sequence, as our first position

% Get all uniprot ids except for the first query sequence
for i=1:numProts
    
    % Need to set up first protein, our query sequence, as our first position
    if i==1
        uniprotx = regexp(proteinMSA(i).Header,'([A-Za-z_-0-9]+)/','tokens'); 
    else
        % Look for uniprot id in header of fasta file
        uniprotx = regexp(proteinMSA(i).Header,'\|([A-Za-z_-0-9]+)\|','tokens');        
    end
    results.uniprotIDs{i} = char(uniprotx{1});
    
    uniprotData = getUniprotInformation('uniprot',results.uniprotIDs{i});
    results.uniprotData{i} = uniprotData;
    if ~isempty(uniprotData)
        results.listOfMRNAs{i} = uniprotData.mrnaEMBL;
        results.mrnaCDS{i} = zeros(length(results.listOfMRNAs{i}),2);
        [results.mrnaSeqs{i},results.mrnaCDS{i},results.mrnaSeqGapped{i},results.mrnaSeqErrors{i},results.mrnaSeqMismatch{i}] = getNucleotideSequencesAndGap(uniprotData.mrnaEMBL,proteinMSA(i).Sequence);
        
        results.listOfGenDNAs{i} = uniprotData.genomicEMBL;
        results.genCDS{i} = zeros(length(results.listOfMRNAs{i}),2);
        [results.genSeqs{i},results.genCDS{i},results.genSeqGapped{i},results.genSeqErrors{i},results.genSeqMismatch{i}] = getNucleotideSequencesAndGap(uniprotData.genomicEMBL,proteinMSA(i).Sequence);
        
    end
end


% Discriminate among sequences to get just one sequence per protein to
% avoid bias
[results.mrnaGapSingle, results.mrnaIDSingle, results.mrnaSeqErrorSingle,results.mrnaChosen] = discriminateAmongTranscripts(results.mrnaSeqGapped, results.listOfMRNAs, results.mrnaSeqErrors);
[results.genGapSingle, results.genIDSingle, results.genSeqErrorSingle,results.genChosen] = discriminateAmongTranscripts(results.genSeqGapped, results.listOfGenDNAs, results.genSeqErrors);

% % Remake aligned sequences into a format useful for writing alignment files
% % and putting into the bioinformatics toolbox
% results.alnNoUTRmRNA = restructureCellSequenceData(results.mrnaGapSingle,results.uniprotIDs,results.mrnaIDSingle);
% results.alnNoUTRGen = restructureCellSequenceData(results.genGapSingle,results.uniprotIDs, results.genIDSingle);
% % 
% % % Look at raw alignments
% [results.mrna5UTR,results.mrna3UTR] = get5and3UTR(results.mrnaSeqs,results.mrnaCDS,results.mrnaChosen);
% [results.gen5UTR,results.gen3UTR] = get5and3UTR(results.genSeqs,results.genCDS,results.genChosen);
% 
% results.alnUTR5mRNA = restructureCellSequenceData(results.mrna5UTR,results.uniprotIDs,results.mrnaIDSingle);
% results.alnUTR3mRNA = restructureCellSequenceData(results.mrna3UTR,results.uniprotIDs,results.mrnaIDSingle);
% results.alnUTR5Gen = restructureCellSequenceData(results.gen5UTR,results.uniprotIDs, results.genIDSingle);
% results.alnUTR3Gen = restructureCellSequenceData(results.gen3UTR,results.uniprotIDs, results.genIDSingle);
% 
% % results.mrnaAli5UTR = multialign(results.alnUTR5mRNA);
% % results.mrnaAli3UTR = multialign(resultsalnUTR3mRNA);
% % results.genAli5UTR = multialign(results.gen5UTR);
% % results.genAli3UTR = multialign(results.gen3UTR);
% 
% fn1 = strcat(results.uniprotIDs{1},'_mrnaCDS.fasta');
% fn2 = strcat(results.uniprotIDs{1},'_genCDS.fasta');
% fastawrite(fn1,results.alnNoUTRmRNA);
% fastawrite(fn2,results.alnNoUTRGen);
% 
save(strcat(results.uniprotIDs{1},'alnRNA'),'results');
end

% function [seq5,seq3] = get5and3UTR(cellSeq,cellCDS)
% 
% N = length(cellSeq);
% seq5 = cell(N,1);
% seq3 = cell(N,1);
% for i=1:N
%     cds = cellCDS{i};
%     seq5{i} = cellSeq{i}(1:cds(1,1)-1);
%     seq3{i} = cellSeq{i}(cds(1,2)+1:end);
% end
% end

% Take in a list of mRNA or genomic DNA sequences, get EMBL sequences, and
% fit them to the protein alignment
function [mrnaSeqs,mrnaCDS,mrnaSeqGapped,mrnaSeqErrors,mrnaSeqMismatch] = getNucleotideSequencesAndGap(listOfMRNAs,proteinMSASeq)

numMRNAs = length(listOfMRNAs);
mrnaSeqs = cell(numMRNAs,1);
mrnaCDS = zeros(numMRNAs,2);
mrnaSeqGapped = cell(numMRNAs,1);
mrnaSeqErrors = zeros(numMRNAs,1);
mrnaSeqMismatch = cell(numMRNAs,1);

% Query EMBL for sequence and cds start and end regions
for j=1:numMRNAs
    
    tempEmbl = getEMBLPlusPlus(listOfMRNAs{j});
    if ~isempty(tempEmbl)
        mrnaSeqs{j} = tempEmbl.Sequence;
        mrnaCDS(j,:) = [tempEmbl.cdsStart tempEmbl.cdsEnd];
        
        [mrnaSeqGapped{j},rnaErrors] = replaceGappedSequenceWithNT(proteinMSASeq,tempEmbl.Sequence(tempEmbl.cdsStart:tempEmbl.cdsEnd));
        mrnaSeqErrors(j) = rnaErrors.error;
        mrnaSeqMismatch{j} = rnaErrors.diffResNt;
    end
end
end


function [gappedRNA,errorCode] = replaceGappedSequenceWithNT(protSeq,mrnaSeq)

% Indices to traverse through our end-result nucleotide answer and our
% given mRNA sequence respectively
idxNt = 1;
idxMRNA = 1;

% Error codes for if the mRNA sequence doesn't match the protein sequence,
% etc.  Error = 0 if fine, 1 if the translated mRNA doesn't match protein
% sequence length, 2 if the mRNA codon doesn't match amino acid at specific
% point shown by diffResNt, 4 if you find gaps in the translated mRNA, and
% 3 if you find gaps in the simple protein sequence
errorCode.error = 0;
errorCode.diffResNt = [];

% Standardize by adding stop codon to protein sequence if not present
if protSeq(end) ~= '*'
    protSeq(end+1) = '*';
end

translMRNA = nt2aa(mrnaSeq,'ACGTOnly','false');
simpleProt = protSeq;
simpleProt((simpleProt == '.') | (simpleProt == '-')) = [];

% Align translated mRNA sequence against queried protein sequence with gaps
% removed.  This will tell where gaps would be in the original query
% sequence, so we can remove them from the transcript before proceeding
[~,newAlign] = nwalign(simpleProt,translMRNA);
if any((newAlign(3,:)=='-') | (newAlign(3,:)=='.'))
    gappedRNA = '';
    errorCode.error = 4;
    return;
end
simpleProtAligned = newAlign(1,:);
z = find((simpleProtAligned == '-') | (simpleProtAligned == '.'));

if ~isempty(z)
    mrnaSeq([(z-1)*3+1 (z-1)*3+2 (z-1)*3+3]) = [];
    errorCode.error = 3;
end


% Our appropriately-gapped RNA sequence
gappedRNA = char(1,3*length(protSeq));

% Check to make sure protein sequence length matches length of given mRNA
if length(simpleProt) ~= length(mrnaSeq)/3
    errorCode.error = 1;
else
    
    % Go through (gapped) protein sequence
    for i=1:length(protSeq)
        aaChar = protSeq(i);
        
        % Insert gaps in mRNA if present in protein
        if (aaChar == '-') || (aaChar == '.')
            
            gappedRNA(idxNt:idxNt+2) = [aaChar aaChar aaChar];
            
            % Reverse-translate amino acid into correct mRNA codon
        else
            gappedRNA(idxNt:idxNt+2) = mrnaSeq(idxMRNA:idxMRNA+2);
            
            
            % Record if nucleotides don't match amino acid:  first column
            % is the index in mRNA transcript, the second is the index in
            % the simplified protein sequence, the third is the aa
            % representation of the protein, and the fourth is the aa of 
            % the mRNA codon
            if ~strcmpi(nt2aa(mrnaSeq(idxMRNA:idxMRNA+2),'ACGTOnly','false'),aaChar)
                
                errorCode.error = 2;
                errorCode.diffResNt = [errorCode.diffResNt; idxMRNA i aa2int(protSeq(i)) aa2int(nt2aa(mrnaSeq(idxMRNA:idxMRNA+2),'ACGTOnly','false'))];
            end
            
            idxMRNA = idxMRNA+3;
        end
        
        % Go to next nucleotide triplet
        idxNt = idxNt+3;
        
    end
end

if ~isrow(gappedRNA)
    gappedRNA = gappedRNA(:)';
end
end

