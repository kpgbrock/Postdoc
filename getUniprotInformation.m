% Takes in a name-value argument pair, either 'filename' or 'uniprot', for
% a filename where each row is a uniprot ID or a single uniprot identifier
% respectively.  Returns a struct containing uniprot id(s), sequence(s),
% length(s), mass(es), secondary structure characteristics and regions 
% (strand, turn, helix), and if given, pdb ids.

function [results] = getUniprotInformation(varargin)

% Find out if uniprot is a cell array or char array
p = inputParser;
addOptional(p,'uniprot','');
addOptional(p,'filename','');

parse(p,varargin{:});

if ~isempty(p.Results.filename)
    
    fileID = fopen(p.Results.filename);
    uniprotList = textscan(fileID,'%s');
    uniprotList = uniprotList{1};
    N = length(uniprotList);
    results.uniprot = cell(N,1);
    results.length = zeros(N,1);
    results.mass = zeros(N,1);
    results.sequence = cell(N,1);
    results.pdb = cell(N,1);
    results.helixLocs = cell(N,1);
    results.turnLocs = cell(N,1);
    results.strandLocs = cell(N,1);
    
    for i=1:N
        [temp] = getUniprotInd(uniprotList{i});
        results.uniprot{i} = temp.uniprot;
        results.mass(i) = temp.mass;
        results.length(i) = temp.length;
        results.sequence{i} = temp.sequence;
        results.pdb{i} = temp.pdb;
        results.helixLocs{i} = temp.helixLocs;
        results.strandLocs{i} = temp.strandLocs;
        results.turnLocs{i} = temp.turnLocs;
    end
    fclose(fileID);
    
elseif ~isempty(p.Results.uniprot)
    results = getUniprotInd(p.Results.uniprot);
else
    disp('Please specify value "uniprot" and a single identifier or "filename" and a filename containing a list of uniprot IDs');
end

end

function [results] = getUniprotInd(uniprot)

results = [];

% Use REST to access XML file for uniprot id
try
    data = webread(strcat('http://www.uniprot.org/uniprot/',uniprot,'.xml'));

% Search through webpage for key information
pattern = '<sequence length="(?<length>\d+).*mass="(?<mass>\d+)".*>\W*(?<sequence>.*)\W*</sequence>';
seqInfo = regexp(data,pattern,'names');
ssPattern1 = '<feature type="';
ssPattern2 = '"[^>]*>\W*<location>\W*<begin position="(?<beginPos>\d+)"/>\W*<end position="(?<endPos>\d+)"/>';
ssStrand = regexp(data,strcat(ssPattern1,'strand',ssPattern2),'names');
ssTurn = regexp(data,strcat(ssPattern1,'turn',ssPattern2),'names');
ssHelix = regexp(data,strcat(ssPattern1,'helix',ssPattern2),'names');
pdbIds = regexp(data,'<dbReference type="PDB" id="(?<pdb>[^"]*)">','names');

% Get EMBL ids for mRNA transcripts
results.mrnaEMBL = regexp(data,'<dbReference type="EMBL" id="([0-9A-Za-z.-]+)">\W*<property type="protein sequence ID" value="[.A-Za-z0-9]+"/>\W*<property type="molecule type" value="mRNA"/>\W*</dbReference>','tokens')';
results.genomicEMBL = regexp(data,'<dbReference type="EMBL" id="([0-9A-Za-z.-]+)">\W*<property type="protein sequence ID" value="[.A-Za-z0-9]+"/>\W*<property type="molecule type" value="Genomic_DNA"/>\W*</dbReference>','tokens')';

% Format output correctly
results.length = str2double(seqInfo.length);
results.mass = str2double(seqInfo.mass);
results.sequence = regexprep(seqInfo.sequence,'\W','');

results.strandLocs = zeros(length(ssStrand),2);
results.turnLocs = zeros(length(ssTurn),2);
results.helixLocs = zeros(length(ssHelix),2);
for i=1:length(ssStrand)
    results.strandLocs(i,1) = str2double(ssStrand(i).beginPos);
    results.strandLocs(i,2) = str2double(ssStrand(i).endPos);
end
for i=1:length(ssTurn)
    results.turnLocs(i,1) = str2double(ssTurn(i).beginPos);
    results.turnLocs(i,2) = str2double(ssTurn(i).endPos);
end
for i=1:length(ssHelix)
    results.helixLocs(i,1) = str2double(ssHelix(i).beginPos);
    results.helixLocs(i,2) = str2double(ssHelix(i).beginPos);
end
    
results.pdb = {pdbIds(:).pdb}';
results.uniprot = uniprot;

catch
    %disp(strcat('Error reading from UniProt', uniprot));
end
    
end