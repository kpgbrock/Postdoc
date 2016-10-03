% Given a directory, see what the current status of different jobs are.
% Takes in a directory name as a string, prints out results on display or
% in optional filename.

function [] = getBuildAli4JobStatus(dirname,outputFile)

output = [];

if ischar(dirname)
    output = getBuildAli4Ind(dirname);
elseif iscell(dirname)
    for i=1:length(dirname)
        output = getBuildAli4Ind(dirname{i},output);
    end
else
    disp('Error:  wrong input type for dirname');
end

if nargin == 2
    writetable(struct2table(output),outputFile);
else
    disp(output);
end

end

function [output] = getBuildAli4Ind(dirname,output)

listfiles = dir(strcat(dirname,filesep,'*.log'));
listfiles = {listfiles.name};
numFiles = length(listfiles);
evals = cell(numFiles,1);
if numFiles == 0
    output.(dirname) = 'No log files found';
end

% Get list of e-values
for i=1:numFiles
    eval = regexp(listfiles{i},'.*_(e[0-9]+).log','tokens');
    evals{i} = strcat(dirname,'_',eval{1});
    
    
    fstring = fileread(strcat(dirname,filesep,listfiles{i}));
    
    % Has it run successfully?
    exCode = strfind(fstring,'Successfully completed.');
    if ~isempty(exCode)
        output.(evals{i}{1}) = 'Successful';
    
    % If not, has it exited with a given exit code?
    else
        exCode = regexp(fstring,'Exited with exit code ([0-9]+)','tokens');
        if ~isempty(exCode)
            output.(evals{i}) = strcat('Exit ',exCode);
    % Is it still running?
        else
            output.(evals{i}) = 'Running';
        end
    
    end
end

end