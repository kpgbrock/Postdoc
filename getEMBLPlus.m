% Wrapper for bioinformatics toolbox getembl function that adds two fields,
% cdsStart and cdsEnd, to the returned data structure.  These correspond to
% the start and end of the protein-coding sequence in mRNA space
% respectively
function [data] = getEMBLPlus(emblID)

try
    
    data = getembl(char(emblID));
    emblFeature = data.Feature;
    
    cdsLine = emblFeature(strmatch('CDS',emblFeature),:);
    getCDSLocs = regexp(cdsLine,'\W+(?<cdsStart>\d+)\.\.(?<cdsEnd>\d+)','names');
    data.cdsStart = str2double(getCDSLocs.cdsStart);
    data.cdsEnd = str2double(getCDSLocs.cdsEnd);
    
catch
    %disp(strcat('Could not access EMBL data',emblID));
    data = [];
    
end

end