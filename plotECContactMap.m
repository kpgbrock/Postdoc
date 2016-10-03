% Takes in as input a list of ECs (variable number of rows, two columns) and the number of
% ECs to plot (N).  Optional arguments are as follows, given as name-value pairs:
%
% 'MarkerSize':  takes in size to make markers, default 72
% 'MarkerType':  character to represent each marker, default 'o'
% 'MarkerEdgeColor':  default is 'flat'
% 'MarkerFaceColor':  default is 'none'
% 'plotLines':  takes in a vector of numbers.  If argument isn't empty, then
%             it plots horizontal and vertical lines corresponding to each
%             number in the vector
% 'plotLinesWidth':  sets width for given input plotted lines, default 3
% 'LineWidth':  changes width of axes, default 6
% 'ecStrengthColors':  takes in a vector of the same length as the list of
%             ECs, and assigns marker colors based on the strength of their
%             couplings
% 'FontSize':  changes font size of axes, default 15
% 'xLabel':  Allows you to set label for x-axis
% 'yLabel':  Ditto for the y-axis
% 'title':  title for figure
% 'removeDiagonal':  Takes in a number N such that you reject ecs that fall
%               within i to i+N


function [] = plotECContactMap(listOfECs, N, varargin)

close all;

% Find out if uniprot is a cell array or char array
p = inputParser;
addOptional(p,'MarkerSize', 72);
addOptional(p,'MarkerType','o');
addOptional(p,'MarkerEdgeColor','flat');
addOptional(p,'MarkerFaceColor','none');
addOptional(p,'plotLines','');
addOptional(p,'plotLinesWidth',2);
addOptional(p,'axesLineWidth',6);
addOptional(p,'ecStrengthColors',[0 0 1]);
addOptional(p,'FontSize',15);
addOptional(p,'xLabel','Residue Number');
addOptional(p,'yLabel','Residue Number');
addOptional(p,'title','');
addOptional(p,'removeDiagonal',0);


parse(p,varargin{:});

removeDiagonalIDs = find(abs(listOfECs(:,1)-listOfECs(:,2)) > p.Results.removeDiagonal);
[r,~] = size(p.Results.ecStrengthColors);
colorinput = p.Results.ecStrengthColors;
if r == length(listOfECs)
    colorinput = p.Results.ecStrengthColors(removeDiagonalIDs);
end
listOfECs = listOfECs(removeDiagonalIDs,:);


z1 = [listOfECs(1:N,1); listOfECs(1:N,2)];
z2 = [listOfECs(1:N,2); listOfECs(1:N,1)];

% If you're putting in ec strengths to make them different colors, make
% sure than you are only taking the top N hits

isColorInput = 0;
%colorinput = p.Results.ecStrengthColors;
if (length(colorinput) == length(listOfECs(:,1)))
    colorinput = [colorinput(1:N); colorinput(1:N)];
    isColorInput = 1;
end

% Plot your EC contact map
h = scatter(z1,z2,p.Results.MarkerSize,colorinput, p.Results.MarkerType,'MarkerFaceColor',p.Results.MarkerFaceColor, 'MarkerEdgeColor',p.Results.MarkerEdgeColor);
set(gca,'LineWidth',p.Results.axesLineWidth,'FontSize',p.Results.FontSize);
xlabel(p.Results.xLabel);
ylabel(p.Results.yLabel);
title(p.Results.title);
if isColorInput
    colorbar;
end

% If lines are expected, plot lines
if ~isempty(p.Results.plotLines)
    hold on;
    numLines = length(p.Results.plotLines);
    ylims = get(gca,'YLim');
    xlims = get(gca,'XLim');
    for i=1:numLines
        plot([p.Results.plotLines(i) p.Results.plotLines(i)],ylims,'k--','LineWidth',p.Results.plotLinesWidth);
        plot(xlims,[p.Results.plotLines(i) p.Results.plotLines(i)],'k--','LineWidth',p.Results.plotLinesWidth);
        uistack(h,'top');
    end
end
    
end
