
function plotThingsOnTopOfGroupedBarPlot(arrayToPlot, avgDim, groupDim, barDim, varargin)

% avgDim is dim to average over (probably the dim that's subjects).
% groupDim is the dim that bar groups are made from (probably ROI)
% barDim is the dim that individual bars within a group are made from
% (probably model)

% varargin can include: 
% a pval array (groupDim x barDim)
% colors for bars (cond x 3)
% colors for individual-subject points (cond x 3)
if numel(varargin) > 0; pvals = varargin{1}; end
if numel(varargin) > 1; barColors = varargin{2}; end
if numel(varargin) > 2; scatterColors = varargin{3}; end

figure; box off; hold on;

b = bar(squeeze(nanmean(arrayToPlot, avgDim))); % Takes mean over the dim we want to mean over (probably 3rd dim if that's 
% axAll = axes('NextPlot','add', 'XTick', [1 2 3], 'XTickLabel', subjIds);


ngroups = size(arrayToPlot, groupDim); % n groups of bars
nbars = size(arrayToPlot, barDim); % n bars per group; probably n models
groupwidth = min(.8, nbars/(nbars + 1.5));
% stdRsAcrossSubj = std(meanRsBySubjAndRoi, [], 1);
% semRsAcrossSubj = stdRsAcrossSubj / sqrt(numel(subjIds));

if numel(varargin) > 1 % if custom colors for bars
    
    for curColor = 1:nbars
        b(curColor).FaceColor = 'flat';
        b(curColor).CData = barColors(curColor,:);
        b(curColor).FaceAlpha = 0.7;
    end
end

highestDotY = max(max(max(arrayToPlot)));
highestAsteriskY = highestDotY .* 1.2; % add buffer before where text would go. 
asteriskIncrement = highestDotY / 20;
asteriskYs = highestAsteriskY - (asteriskIncrement .* [0 1 2]);
topYlim = asteriskYs(1) + asteriskIncrement;

% This part is specific to situations w lower signal for second 3 bar
% groups (ROIs) than the rest
highestSceneY = max(max(max(arrayToPlot(3:6,:,:))));
highestAsteriskSceneY = highestSceneY + (highestDotY - asteriskYs(1)); % Make even spacing between evc and scene-area points and asterisks
asteriskSceneYs = highestAsteriskSceneY - (asteriskIncrement .* [0 1 2]);

asteriskFontSize = 14; 
asteriskFontWeight = 'bold';
nsFontSize = 8;
nsFontWeight = 'regular';

scatterSize = 15;

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange

% withinSubjSems = semWithin(meanRsBySubjAndRoi);
for i = 1:nbars
    % Calculate center of each bar
    % (Will be able to tell if sem is wrong, bc order should be wrong for
    % both)
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    if i==1; firstX = x; end
    
    jitterWidth = groupwidth/12;
    jitterRange = linspace(-0.5, 0.5, 8);
    jitterVals = jitterRange .* jitterWidth; 
    
    if barDim == 2 % e.g., if bars within a group are models.
        for iGroup = 1:ngroups
            individSubjYsThisBar = arrayToPlot(iGroup,i,:);
            individSubjXsThisBar = ...
                (ones(size(arrayToPlot,3),1) .* x(iGroup)) + jitterVals';
            sHandle = scatter(individSubjXsThisBar, individSubjYsThisBar, ...
                scatterSize, 'o');
            if numel(varargin) > 1
                sHandle.MarkerEdgeColor = scatterColors(i,:);
            end
            
            if exist('pvals','var')
                if mean(individSubjYsThisBar) < 0
                    starText = '';
                elseif pvals(iGroup,i) < 0.0001
                    starText = '****';
                    curFontSize = asteriskFontSize;
                    curFontWeight = asteriskFontWeight;
                elseif pvals(iGroup,i) < 0.001
                    starText = '***';
                    curFontSize = asteriskFontSize;
                    curFontWeight = asteriskFontWeight;
                elseif pvals(iGroup,i) <= 0.01
                    starText = '**';
                    curFontSize = asteriskFontSize;
                    curFontWeight = asteriskFontWeight;
                elseif pvals(iGroup,i) <= 0.05
                    starText = '*';
                    curFontSize = asteriskFontSize;
                    curFontWeight = asteriskFontWeight;
                else
                    % starText = sprintf('p = %0.2f', pvals(iGroup,i));
                    starText = '';
                    curFontSize = nsFontSize;
                    curFontWeight = nsFontWeight;
                end
                if exist('scatterColors','var')
                    starColor = scatterColors(i,:);
                elseif exist('barColors','var')
                    starColor = barColors(i,:);
                else
                    starColor = [0 0 0];
                end
                if iGroup < 4 || iGroup > 6
                    asteriskYsToUse = asteriskYs;
                else
                    asteriskYsToUse = asteriskSceneYs;
                end
                text(firstX(iGroup), asteriskYsToUse(i), starText, ...
                    'Color',starColor, 'FontSize', curFontSize, 'FontWeight', 'bold');
            end
            
        end
    end
    %{
    errorbar(x, mean(meanRsBySubjAndRoi(:,i),1), ...
        withinSubjSems(i), 'Color',[0.2 0.2 0.2], ...
        'CapSize',3,'linestyle', 'none'); % 'Color',pairedColorMap(i,:),
    %}
    
    
end

ylim([-0.01 topYlim]);

%{
% Don't need to reorder these bc they're just the colors
for iRoi = 1:numel(roiNamesToUse)
    b(iRoi).EdgeColor = colorList(iRoi,:);
    b(iRoi).FaceColor = colorList(iRoi,:);
    b(iRoi).FaceAlpha = .7;
end
%}

% xtickangle(-45);

% axXYOnly = axes('NextPlot','add', 'XColor',[0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8],...
    % 'YTickLabel', {}, 'XTickLabel', {}, 'XTick', []);

% xlim([0.5 1.5]);
% ax.XRuler.Axle.LineWidth = 2;
% axXYOnly.YRuler.Axle.LineWidth = 3;


end