function [data] = getEMBLPlusPlus(emblID)

    %pause(0.1);
    
    % Force your input to be in char form
    emblID = char(emblID);
    emblPrefix = 'http://www.ebi.ac.uk/ena/data/view/';
    emblSuffixXML = '&display=xml';
    emblSuffixFasta = '&display=fasta';
    
    emblURLXML = strcat(emblPrefix,emblID,emblSuffixXML);
    emblURLFasta = strcat(emblPrefix,emblID,emblSuffixFasta);
    
    data.cdsStart = 0;
    data.cdsEnd = 0;
    
   try
        % Get sequence
        data.Sequence = '';
        xmldata = webread(emblURLXML);
        
        rawSeq = regexp(xmldata,'<sequence>(.*)</sequence>','tokens');
        
        if isempty(rawSeq)
            disp('isempty')
            
            seqdata = webread(emblURLFasta);
            
            seqdata = fastaread(seqdata);
            data.Sequence = seqdata.Sequence;
        else
            
            data.Sequence = rawSeq{1}{1};
            data.Sequence(isspace(data.Sequence)) = [];
        end
       
   catch
        disp(strcat('Error reading sequence for EMBL: ',emblID));
        data = [];
        return;
   end

    cdsLine = regexp(xmldata,'<feature name="CDS" location="(\d+)..(\d+)">','tokens');
    if ~isempty(cdsLine)
        data.cdsStart = str2double(cdsLine{1}{1});
        data.cdsEnd = str2double(cdsLine{1}{2});
    else
        disp(strcat('Coding start and end ambiguous: ',emblID));
        data.cdsStart = 1;
        data.cdsEnd = length(data.Sequence);
    end
    
end