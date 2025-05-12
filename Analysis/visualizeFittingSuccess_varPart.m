% 2022.02.07 AS made more generalizable to be able to switch out feat types
% smoothly, as well as able to read in thresholded VOIs.

%% Setup
clear; close all;

analysisDir = pwd; addpath(genpath('helpers')); 
filesepinds = find(analysisDir==filesep);
mainDir = analysisDir(1:filesepinds(end)-1);
fmriDataDir = fullfile(mainDir,'DataFmri');
intermedDir = fullfile(analysisDir,'ridgeFigs','varPartIntermedArrays');
saveIntermedFiles = false; 
useSavedIntermedFiles = true; % See above comment
% doFittingByTask = false;

subjIds = {'CI','CO','CK','CR','BX','CS','CG','CH'}; % ,'CH'};
subjNumsStr = {'S01','S02','S03','S04','S05','S06','S07','S08'}; % ,'S08'};
subjsToDoNow = 1:8;

% thresholdedVois = 1; % This is the official option bc means the voxels passed a localizer threshold in independent localizer data.
plotThresholdVoiDetails = false;
%{
gssParcSceneAreas = 0;
if thresholdedVois && gssParcSceneAreas
    error('Oh no! thresholdedVois and gssParcSceneAreas are mutually exclusive options!');
end
if gssParcSceneAreas == 1
    % which areas do we want to load voxels from?
    roiNames = {'V1','V2','V3', ...
        'OPA_gssParc', 'PPA_gssParc', 'RSC_gssParc'}; % ,'V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    % 'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA'};
    sceneText = '_gssSceneParcs';
elseif thresholdedVois == 1
    %}
    roiNames = {'V1','V2','V3', ...
        'OPA', 'PPA', 'RSC'}; % ,'V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    % 'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA'};
    % sceneText = '_thresh';
    %{
elseif gssParcSceneAreas == 0
    % which areas do we want to load voxels from?
    roiNames = {'V1','V2','V3','OPA_manual','PPA_manual','RSC_manual'};
    sceneText = '_manualSceneAreas';
end
    %}
roisToDoNow = 1:6; % Start with just EVC

% newRidgeCode = 1;
% fromPython = 0;
% doublePrec = 1;
% noDepthArtifactIms = 1;
% veryFewLambdas = 1;
lookAtPerms = true; % 1;
collapsedHemis = 1;
% intercept = 0; % This should be 0!

% artifactPrefix = '';
% lambdaText = '_fewerLambdas';

%{
if lookAtPerms
    permText = '_permTests';
else
    permText = '';
end
%}
%{
if intercept==0
    intText = '_noIntercept';
else
    intText = '';
end
%}

featLabels = {'spat3d','nonSpat3d','gabor'}; % also a pix model, but probably not involved in var part analysis.
featTypesA = [1, 2];
featTypesB = [3, 3];

for iFeatCombo = 1:numel(featTypesA)
    featTypeA = featTypesA(iFeatCombo);
    featTypeB = featTypesB(iFeatCombo);

    featPrefixA = ['_' featLabels{featTypeA} '_']; % e.g., '_nonSpat3d_'
    featPrefixB = ['_' featLabels{featTypeB} '_']; % e.g., '_gabor_'

% ridgeOutDir = fullfile(analysisDir, 'ridgeOutputs', ...
    % ['varPart' featPrefixA featPrefixB artifactPrefix lambdaText permText]);
ridgeOutDir = fullfile(analysisDir, 'ridgeOutputs', ...
        ['varPart' featPrefixA featLabels{featTypeB}]);

hemiText = '_collapsedHemispheres';
roiNamesToUse = roiNames;
newRoiOrders = 1:6; % [1 6 2 7 3 8 4 9 5 10];
colorIndsOldOrders = 1:6;
%{
pairedColorMap = [ ...
    .8 .4 .6; ... % V1
    .9 .3 .4; ... % V2
    .9 .35 .3; ... % V3
    .6 .7 .4; ... % OPA
    .3 .65 .7; ... % PPA
    .4 .35 .7]; % RSC
%}
%{
if newRidgeCode==1
    versionText = '_vectorized';
elseif newRidgeCode==0
    versionText = '';
end
if doublePrec
    precText = '_double_';
elseif doublePrec == 0
    precText = '_single_';
end
%}

% ------------------------------------
% Time window options 
% ------------------------------------

% Avg. sig. analysis (meaning averaged over one time window only; not mutually exclusive w analysis over time, below)
doAvgSignalAnalysis = 1; 
longerTrWindow = 0;
% shorterTrWindow = 1; % official analysis window afaik
% trWindow2 = 0; % I believe this is the option attempting to roughly match 
% gallant-lab and BOLD5000 stuff, although won't correspond perfectly since 
% this experiment uses a fast event-related design.
%{
% ... And built-in error messages to make sure options set correctly
if (longerTrWindow + shorterTrWindow + trWindow2) > 1
    error('Oh no! You can''t run script for multiple TR windows at once!');
elseif ~doAvgSignalAnalysis && ((longerTrWindow + shorterTrWindow + trWindow2) > 0)
    error('Oh no! You must have doAvgSignalAnalysis turned on to use sub-options.');
end

% Output text and TRs to grab if doing TR-by-TR analysis
if longerTrWindow % Controls which main signal file gets read as well as the ridge output file name
    avgTRs_targ = [6,10];
    trWindowText = '4-8sec';
elseif shorterTrWindow
    avgTRs_targ = [4,6];
    %}
    trWindowText = '2.4-4.8sec';
    
    %{
elseif trWindow2
    avgTRs_targ = [7,9];
    trWindowText = '4.8-7.2sec';
else
    % avgTRs_targ = [6,8]; % Will end up including 4-6.4 seconds after stim; includes times from future trials, as in Lescroart & Gallant (2019).
    % trWindowText = '4-6ishSec';
    trWindowText = '2.4-4.8sec'; % If doing thresholded vois (main) and not avg sig analysis, this should be the current standard.
    end
%}

% Note: consider making these into 3d arrays later instead of separate
% arrays.
rSqsABySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse));
rSqsBBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse));
rSqsConcatBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse));
rSqsAUBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse));
rSqsBUBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse));
rSqsSharedBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse));
if lookAtPerms
    nullRSqsBySubjAndRoiAndVarType = cell(numel(subjIds), numel(roiNamesToUse), 3); % 3 var types
end

sceneLocTvalsBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse)); % scenes > objects
evcLocTvalsBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse)); % scrambled > baseline

%{
meanRSqsAUBySubjAndRoi = nan(numel(subjIds),numel(roiNamesToUse));
meanRSqsBUBySubjAndRoi = nan(size(meanRSqsAUBySubjAndRoi));
meanRSqsSharedBySubjAndRoi = nan(size(meanRSqsAUBySubjAndRoi));
% medRsBySubjAndRoi = nan(size(meanRsBySubjAndRoi));
semRSqsAUBySubjAndRoi = nan(size(meanRSqsAUBySubjAndRoi));
semRSqsBUBySubjAndRoi = nan(size(meanRSqsAUBySubjAndRoi));
semRSqsSharedBySubjAndRoi = nan(size(meanRSqsAUBySubjAndRoi));

prctile975AUBySubjAndRoi = nan(size(meanRSqsAUBySubjAndRoi));
prctile975BUBySubjAndRoi = nan(size(meanRSqsAUBySubjAndRoi));
prctile975SharedBySubjAndRoi = nan(size(meanRSqsAUBySubjAndRoi));
%}

%% Make/fill summaries; plot hists for each ROI/subj
%  NOT REORDERING BY ROI YET IN THIS SECTION (actually for taskAttScene,
%  don't have to reorder at all; this was left over from BOLD5000 data.)

% voxDistsFig = figure; hold on; % For all distributions of fitting success across voxels
for iSubj = 1:numel(subjIds)
    curSubj = subjIds{iSubj};
    curSubjNum = subjNumsStr{iSubj};
    fprintf('Starting %s \n',curSubj);
    % figure; 

    for iRoi = 1:numel(roiNamesToUse)
        curRoi = roiNamesToUse{iRoi};
        % if ~fromPython
            % curFName = fullfile( ...
                % ridgeOutDir, ['fittedModelsConcat' featPrefixA featPrefixB curSubj '_' curRoi sceneText '_double__vectorized_collapsedHemispheres.mat']);
            curFName = fullfile( ...
                ridgeOutDir, ['fittedModelsConcat' featPrefixA featPrefixB curSubj '_' curRoi '_collapsedHemispheres.mat']);
            
            load(curFName, 'curPredRsA','curPredRsB','curPredRsConcat'); %'curPredRs*');
            if lookAtPerms
                load(curFName, 'curPredRsA_perms', 'curPredRsB_perms', 'curPredRsConcat_perms');
            end
            
            
            % if thresholdedVois
                curLocTvalFName = fullfile(fmriDataDir, ...
                    ['MainTaskSignalByTrial_' trWindowText ... sceneText ...
                    '_' curSubjNum, '.mat']);
                load(curLocTvalFName, 'localizerInfo');
            % end
            
            %{
            if (featType==1) || (featType==2)
                if newRidgeCode==1
                    load(curFName, 'meanRsByVox');
                elseif newRidgeCode == 0
                    load(curFName, 'meanRsByVox_threeD');
                    meanRsByVox = meanRsByVox_threeD;
                end
            elseif featType == 3
                load(curFName, 'meanRsByVox');
            end
            
            if lookAtPerms
               load(curFName,'meanRsByVox_perms'); 
            end
            %}
            rSqsABySubjAndRoi{iSubj,iRoi} = ...
                sign(curPredRsA) .* (curPredRsA .^ 2);
            rSqsBBySubjAndRoi{iSubj,iRoi} = ...
                sign(curPredRsB) .* (curPredRsB .^2);
            rSqsConcatBySubjAndRoi{iSubj,iRoi} = ...
                sign(curPredRsConcat) .* (curPredRsConcat .^2);
            
            
            rSqsAUBySubjAndRoi{iSubj,iRoi} = mean(...
                rSqsConcatBySubjAndRoi{iSubj,iRoi} - ...
                rSqsBBySubjAndRoi{iSubj,iRoi}, 2); % Taking mean across cv folds
            rSqsBUBySubjAndRoi{iSubj,iRoi} = mean(...
                rSqsConcatBySubjAndRoi{iSubj,iRoi} - ...
                rSqsABySubjAndRoi{iSubj,iRoi}, 2);
            rSqsSharedBySubjAndRoi{iSubj,iRoi} = mean( ...
                rSqsABySubjAndRoi{iSubj,iRoi} + ...
                rSqsBBySubjAndRoi{iSubj,iRoi} - ...
                rSqsConcatBySubjAndRoi{iSubj,iRoi}, 2);
            
            %{
            % rsBySubjAndRoi{iSubj,iRoi} = meanRsByVox;
            if ~plotThresholdVoiDetails % This means across all vox; doesn't threshold or rank-order based on independent-localizer activation.
                meanRSqsAUBySubjAndRoi(iSubj,iRoi) = mean(rSqsAUBySubjAndRoi{iSubj,iRoi}); % Mean r-squareds across vox
                meanRSqsBUBySubjAndRoi(iSubj,iRoi) = mean(rSqsBUBySubjAndRoi{iSubj,iRoi});% medRsBySubjAndRoi(iSubj,iRoi) = median(meanRsByVox);
                meanRSqsSharedBySubjAndRoi(iSubj,iRoi) = mean(rSqsSharedBySubjAndRoi{iSubj,iRoi});
            end
            %}
            
            if lookAtPerms
                curARSqPermArray = sign(curPredRsA_perms) .* ...
                    (curPredRsA_perms .^2);
                curBRSqPermArray = sign(curPredRsB_perms) .* ...
                    (curPredRsB_perms .^2);
                curConcatRSqPermArray = sign(curPredRsConcat_perms) .* ...
                    (curPredRsConcat_perms .^2);
                
                curAURSqPermArray = squeeze(mean( ...
                    (curConcatRSqPermArray - ...
                    curBRSqPermArray), 2));
                curBURSqPermArray = squeeze(mean( ...
                    (curConcatRSqPermArray - ...
                    curARSqPermArray), 2));
                curSharedRSqPermArray = squeeze(mean( ...
                    (curARSqPermArray + curBRSqPermArray - ...
                    curConcatRSqPermArray),2));
                
                % For filling array: all are ordered so BU comes first, 
                % then shared, then AU. (Can verify in section
                % concatenating arrays for plotting)
                nullRSqsBySubjAndRoiAndVarType{iSubj,iRoi,1} = ...
                    curBURSqPermArray;
                nullRSqsBySubjAndRoiAndVarType{iSubj,iRoi,2} = ...
                    curSharedRSqPermArray;
                nullRSqsBySubjAndRoiAndVarType{iSubj,iRoi,3} = ...
                    curAURSqPermArray;
                    
            end
            
            % Info about localizer t-vals (for voxel selection)
            sceneLocTvalsBySubjAndRoi{iSubj,iRoi} = localizerInfo(iRoi).sceneTvals; % scenes > objects
            evcLocTvalsBySubjAndRoi{iSubj,iRoi} = localizerInfo(iRoi).evcTvals; % scrambled > baseline

            %{
            semRSqsAUBySubjAndRoi(iSubj,iRoi) = std(rSqsAUBySubjAndRoi{iSubj,iRoi}) / ...
                sqrt(length(rSqsAUBySubjAndRoi{iSubj,iRoi}));
            semRSqsBUBySubjAndRoi(iSubj,iRoi) = std(rSqsBUBySubjAndRoi{iSubj,iRoi}) / ...
                sqrt(length(rSqsBUBySubjAndRoi{iSubj,iRoi}));
            semRSqsSharedBySubjAndRoi(iSubj,iRoi) = std(rSqsSharedBySubjAndRoi{iSubj,iRoi}) / ...
                sqrt(length(rSqsSharedBySubjAndRoi{iSubj,iRoi}));
            %}
        % else % if fromPython
            % error('Oh no! No python version of variance partitioning is ready.');
            %{
            curFName = fullfile( ...
                featDir, ['corrs_' curSubj '_' featLabels{featType} '.h5']);
            
            curCorrs = mean(h5read(curFName,['/corrs' curRoi]),1); % Got transposed somehow between python and here ...? but first dim is cv folds.
            rsBySubjAndRoi{iSubj,iRoi} = curCorrs;
            meanRsBySubjAndRoi(iSubj,iRoi) = mean(curCorrs);
            medRsBySubjAndRoi(iSubj,iRoi) = median(curCorrs);
            semRsBySubjAndRoi(iSubj,iRoi) = std(curCorrs) / sqrt(length(curCorrs));
            %}
        % end
        % set(0, 'currentFigure', voxDistsFig);
        % curPlot = ((iSubj-1)*numel(roiNamesToUse)) + iRoi;
        % subplot(numel(subjIds),numel(roiNamesToUse), curPlot); hold on; box off;
        %{
        if ~fromPython
            % histogram(meanRsByVox,'FaceColor',[0.1 0.1 0.1]);
            % histogram(rSqsAUBySubjAndRoi{iSubj,iRoi},'FaceColor',[0.1 0.1 0.1]); % title(sprintf('non-spat3d, %s, %s
            % histogram(rSqsBUBySubjAndRoi{iSubj,iRoi},'FaceColor',[0.1 0.1 0.1]);
            % histogram(rSqsSharedBySubjAndRoi{iSubj,iRoi},'FaceColor',[0.1 0.1 0.1]);
            
        elseif fromPython
            histogram(curCorrs,'FaceColor',[0.1 0.1 0.1]);
        end
        %}
        % zPlot = plot([0 0],[0 300],'Color',[0 0 0], 'LineWidth',3);
        % zPlot.Color(4) = .2;
        
        % if fromPython
            % title(sprintf('meanRs from Python,\n %s,\n %s', curSubj, curRoi));
            % clear curCorrs
        % else
            % title(sprintf('var partitioning RSqs,\n %s,\n %s', curSubj, curRoi));
            % ylim([0 300]);
            % xlim([-0.05 0.05]);
            
            %{
            % COMMENTING THIS OUT TILL UPDATING LABELS TO NOT BE HARD-CODED 
            if lookAtPerms
                plotNullDistAndStats(curAURSqPermArray, ...
                    rSqsAUBySubjAndRoi{iSubj,iRoi}, ...
                    pairedColorMap,colorIndsOldOrders,iRoi,curRoi,curSubj, ...
                    'quadrant-based 3d, unique');
                saveas(gcf,fullfile('ridgeFigs','varPart','statsHists','quadrant-based3dUnique', ...
                    sprintf('statsHist_quadrant-based3dUnique_%s_%s.pdf',curSubj,curRoi)));
                plotNullDistAndStats(curBURSqPermArray, ...
                    rSqsBUBySubjAndRoi{iSubj,iRoi}, ...
                    pairedColorMap,colorIndsOldOrders,iRoi,curRoi,curSubj, ...
                    'global 3d, unique');
                saveas(gcf,fullfile('ridgeFigs','varPart','statsHists','global3dUnique', ...
                    sprintf('statsHist_global3dUnique_%s_%s.pdf',curSubj,curRoi)));
                plotNullDistAndStats(curSharedRSqPermArray, ...
                    rSqsSharedBySubjAndRoi{iSubj,iRoi}, ...
                    pairedColorMap,colorIndsOldOrders,iRoi,curRoi,curSubj, ...
                    'shared');
                saveas(gcf,fullfile('ridgeFigs','varPart','statsHists','shared', ...
                    sprintf('statsHist_shared_%s_%s.pdf',curSubj,curRoi)));
            
            end
            %}
            
            %{
            plotNullDistAndStats(squeeze(mean(curPredRsA_perms,2)), ...
                rSqsAUBySubjAndRoi{iSubj,iRoi}, ...
                pairedColorMap,colorIndsOldOrders,iRoi,curRoi,curSubj, ...
                '3d surface, unique');
                %}
            %{
            
            %}
            % clear rSqs*
        % end
    end
end

            
%% Plot summaries
%{
figure; 
if noDepthArtifactIms==1
    title('Prediction perf. (bad ims excluded)');
else
    title('Prediction perf. (dumber cv loops, no exclusions)');
end
box off; hold on;

b = bar(meanRsBySubjAndRoi(:,newRoiOrders));
% axAll = axes('NextPlot','add', 'XTick', [1 2 3], 'XTickLabel', subjIds);

ngroups = size(meanRsBySubjAndRoi, 1);
nbars = size(meanRsBySubjAndRoi, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    % (Will be able to tell if sem is wrong, bc order should be wrong for
    % both)
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, meanRsBySubjAndRoi(:,newRoiOrders(i)), ...
        semRsBySubjAndRoi(:,newRoiOrders(i)), 'Color',[0.2 0.2 0.2], ...
        'CapSize',3,'linestyle', 'none'); % 'Color',pairedColorMap(i,:),
end

% Don't need to reorder these bc they're just the colors
for iRoi = 1:numel(roiNamesToUse)
    b(iRoi).EdgeColor = pairedColorMap(iRoi,:);
    b(iRoi).FaceColor = pairedColorMap(iRoi,:);
    b(iRoi).FaceAlpha = .7;
end

legend(roiNamesToUse(newRoiOrders));
ylabel('Mean of corrs w held-out data');
xlabel('Subj');
xticklabels(subjIds);

% axXYOnly = axes('NextPlot','add', 'XColor',[0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8],...
    % 'YTickLabel', {}, 'XTickLabel', {}, 'XTick', []);

xticks([1 2 3]); xticklabels(subjIds);
xlim([0.5 4.5]);
%}
% ax.XRuler.Axle.LineWidth = 2;
% axXYOnly.YRuler.Axle.LineWidth = 3;

%% Concatenate different variance types into an array together for plotting.
% B-unique is always going to go first! Then shared, then A-unique
% These are not subsetted by voxel; leaving here bc subsetting happens
% later.
%{
subjMeanVarPart = [ ...
    mean(meanRSqsBUBySubjAndRoi)', ... % Across-subj mean of across-vox means
    mean(meanRSqsSharedBySubjAndRoi)',...
    mean(meanRSqsAUBySubjAndRoi)'];

subjMeanVarPartSems = [ ...
    semWithin(meanRSqsBUBySubjAndRoi)', ... % within-subj sem of across-vox means
    semWithin(meanRSqsSharedBySubjAndRoi)',...
    semWithin(meanRSqsAUBySubjAndRoi)'];

% Note that the order goes B, A, but labels are right below so you can
% verify they correspond correctly.
arraysOfVarsInOrder = {...
    meanRSqsBUBySubjAndRoi, meanRSqsSharedBySubjAndRoi, meanRSqsAUBySubjAndRoi};
varTypeLabels = ...
    {[featLabels{featTypeB} ' unique'], ...
    'shared', ...
    [featLabels{featTypeA} ' unique']};

% Also manually ordered as BU, shared, AU
rSqsBySubjAndRoiAndVarType = cell(numel(subjIds), numel(roiNamesToUse));
for iSubj = 1:numel(subjIds)
    for iRoi = 1:numel(roiNamesToUse)
        rSqsBySubjAndRoiAndVarType{iSubj,iRoi}(:,1) = ...
            rSqsBUBySubjAndRoi{iSubj,iRoi};
        rSqsBySubjAndRoiAndVarType{iSubj,iRoi}(:,2) = ...
            rSqsSharedBySubjAndRoi{iSubj,iRoi};
        rSqsBySubjAndRoiAndVarType{iSubj,iRoi}(:,3) = ...
            rSqsAUBySubjAndRoi{iSubj,iRoi};
    end
end
%}
        

%% Plot summaries, subj-mean plots! (Concatenate, make grouped!)
% (NOTE: THIS IS ONLY IF NOT SELECTING VOX ABOVE CUTOFF!)
% (so only preliminary stages!!!)
%{
if ~plotThresholdVoiDetails
    figure;
    %{
if noDepthArtifactIms==1
    title('Across-subj means (bad ims excluded)');
else
    title('Across-subj means (dumber cv loops, no exclusions)');
end
    %}
    box off; hold on;
    
    b = bar(subjMeanVarPart(newRoiOrders,:));
    xticklabels(roiNamesToUse(newRoiOrders));
    
    ngroups = 6;
    nbars = size(subjMeanVarPart, 2);
    groupwidth = min(.8, nbars/(nbars + 1.5));
    % stdRsAcrossSubj = std(meanRsBySubjAndRoi, [], 1);
    % semRsAcrossSubj = stdRsAcrossSubj / sqrt(numel(subjIds));
    
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    
    % withinSubjSems = semWithin(meanRsBySubjAndRoi);
    
    for i = 1:nbars
        % Calculate center of each bar
        % (Will be able to tell if sem is wrong, bc order should be wrong for
        % both)
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, subjMeanVarPart(newRoiOrders,i), ...
            subjMeanVarPartSems(newRoiOrders,i), 'Color',[0.2 0.2 0.2], ...
            'CapSize',3,'linestyle', 'none'); % 'Color',pairedColorMap(i,:),
        
        for iSubj = 1:numel(subjIds)
            scatter(x, arraysOfVarsInOrder{i}(iSubj,:), 16, [0.1 0.1 0.1], ...
                '.')
        end
    end
    
    
    % Don't need to reorder these bc they're just the colors
    grays = [ ...
        0.7 0.7 0.7;
        0.55 0.55 0.55;
        0.42 0.42 0.42];
    
    for iVarType = 1:nbars
        % b(iRoi).EdgeColor = grays(iRoi,:);
        b(iVarType).FaceColor = grays(iVarType,:);
        % b(iRoi).FaceAlpha = .7;
    end
    
    % legend(roiNamesToUse(newRoiOrders));
    % ylabel('Mean of corrs w held-out data');
    xlabel('ROI');
    xticklabels(roiNames(newRoiOrders));
    set(gca,'xticklabel',{[]})
    ylabel('Propor variance explained, r-squared');
    legend(b,varTypeLabels);
    if thresholdedVois
        title('Includes all vox above liberal localizer threshold');
    end
    % ylim([0 0.1]);
    % ylim([0 0.09]);
    % xtickangle(-45);
    
    % axXYOnly = axes('NextPlot','add', 'XColor',[0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8],...
    % 'YTickLabel', {}, 'XTickLabel', {}, 'XTick', []);
    
    % xlim([0.5 1.5]);
    % ax.XRuler.Axle.LineWidth = 2;
    % axXYOnly.YRuler.Axle.LineWidth = 3;
end
%}

%% Voxels in order of loc activity (individual subjects); 
% This section plots individual subjects/models at a time (saving info
% in intermediate files to be able to plot all of them together in a summary)
% BUT ALSO: saves non-null and null arrays to be plotted later in subect
% summary after selecting the correct voxels to use based on localizer activation.

if plotThresholdVoiDetails
    
    % This is for the plot that only uses the top 200 (or other) voxels.
    rSqToPlotEachRoiVarTypeAndSubj = ...
        nan(numel(roiNamesToUse), numel(varTypeLabels), numel(subjIds));
    nullRSqArrayFinalVoxEachRoiVarTypeAndSubj = ...
        cell(numel(roiNamesToUse), numel(varTypeLabels), numel(subjIds));
    % This is a corresponding array to say how many voxels used
    nVoxBySubjAndRoi = nan(numel(subjIds), numel(roiNamesToUse));
    roiIsUsableBySubjAndRoi = nan(size(nVoxBySubjAndRoi));
    
    for iSubj = 1:numel(subjIds)
        % scatterFig = figure; % Temporarily not working!
        voxSelectionPlot = figure;
        plotCtr = 0;
        curSubjNum = subjNumsStr{iSubj};
        curSubj = subjIds{iSubj};
        for iRoi = 1:numel(roiNamesToUse)
            curRoi = roiNamesToUse{iRoi};
            plotCtr = plotCtr + 1;
            
            [sortedEvcTvals, highestEvcTvalInds, ...
                sortedSceneTvals, highestSceneTvalInds] = ...
                getSortedTvalsFromArrays( ...
                evcLocTvalsBySubjAndRoi{iSubj,iRoi}, ...
                sceneLocTvalsBySubjAndRoi{iSubj,iRoi});
            
            % Scatter figure ( jk have to re-fill that in bc accidentally
            % deleted.) It was prediction performance plotted against 
            % localizer t-val.
            % -------------------------------------------------------
            
            
            % Figure showing how/whether effects vary depending on number
            % of voxels included
            % -------------------------------------------------------
            
            set(0,'CurrentFigure',voxSelectionPlot);
            if iRoi < 4 % Use evc tvals as inputs for EVC ROIs
                [curArrayOfAvgsByVarType, ~] = plotVoxSelectionStackedPlotWithinLoop(...
                    rSqsBySubjAndRoiAndVarType, ...
                    sortedEvcTvals, highestEvcTvalInds, ...
                    varTypeLabels, iRoi, iSubj, curRoi, curSubjNum, plotCtr);
            else % Use scene tvals as inputs for scene-selective ROIs
                [curArrayOfAvgsByVarType, ~] = plotVoxSelectionStackedPlotWithinLoop(...
                    rSqsBySubjAndRoiAndVarType, ...
                    sortedSceneTvals, highestSceneTvalInds, ...
                    varTypeLabels, iRoi, iSubj, curRoi, curSubjNum, plotCtr);
            end
            
            % Note: script is currently set up to run this only for the one
            % feature type being looked at rn. You can change the feat type
            % at the top and run again if you want all of these.
            %{
            curIntermedFName = ...
                    fullfile(intermedDir, ...
                    ['rSqsSortedByVoxLocActivation' curSubj '_' curRoi precText ...
                    versionText hemiText intText '_' trWindowText '_' ...
                    sceneText '.mat']);
                %}
                curIntermedFName = ...
                    fullfile(intermedDir, ...
                    ['rSqsSortedByVoxLocActivation' curSubj '_' curRoi  ...
                    hemiText '_' trWindowText '_.mat']);
            if saveIntermedFiles
                save(curIntermedFName, ...
                    'curSubj','curSubjNum','iSubj','curArrayOfAvgsByVarType', ...
                    'varTypeLabels','iRoi','curRoi', ...
                    'sortedEvcTvals', 'highestEvcTvalInds', ...
                    'sortedSceneTvals', 'highestSceneTvalInds');
                clear curIntermedFName
            elseif useSavedIntermedFiles
                load(curIntermedFName, ...
                    'curSubj','curSubjNum','iSubj','curArrayOfAvgsByVarType', ...
                    'varTypeLabels','iRoi','curRoi', ...
                    'sortedEvcTvals', 'highestEvcTvalInds', ...
                    'sortedSceneTvals', 'highestSceneTvalInds');
                clear curIntermedFName
            end
            
            % Save info about how many vox used here, as well as the
            % r-sq-value corresponding to that number of vox.
            % (Redundant to do some of this for each model, but the
            % info's only read in from files saved for each model, so
            % just doing it redundantly)
            maxVoxToUse = 200; % Max vox to average over
            voxCutoff = 75; % Min vox to be included in analysis.
            curNVoxUsed = min(maxVoxToUse, size(curArrayOfAvgsByVarType,1)); % Either use number of vox (if less than 200) or use 200 (cutoff)
            nVoxBySubjAndRoi(iSubj,iRoi) = curNVoxUsed;
            roiIsUsableBySubjAndRoi(iSubj,iRoi) = curNVoxUsed >= voxCutoff;
            if iRoi < 4 % To index into null array. Kind of frankenstein-y here but if I streamlined code, this could also be used to index into non-null arrays
                curFinalVoxInds = highestEvcTvalInds(1:curNVoxUsed);
            else
                curFinalVoxInds = highestSceneTvalInds(1:curNVoxUsed);
            end
            
            for iVarType = 1:numel(varTypeLabels) % This loop might just be left over from when I already had a loop for a different reason.
                if roiIsUsableBySubjAndRoi(iSubj,iRoi) % Fill in if usable; otherwise leave NaN
                    rSqToPlotEachRoiVarTypeAndSubj(iRoi,iVarType,iSubj) = ... % If I were paring this down, I'd use the array that's just by voxel, not top-N voxels
                        curArrayOfAvgsByVarType(curNVoxUsed, iVarType); % curNVoxUsed is the number of voxels to average over, and curArrayOfAvgs is ordered by how many vox to avg over.
                end
                
                if lookAtPerms % nullRSqsBySubjAndRoiAndVarType
                    nullRSqArrayFinalVoxEachRoiVarTypeAndSubj{iRoi,iVarType,iSubj} = ...
                        nullRSqsBySubjAndRoiAndVarType{iSubj,iRoi,iVarType}(curFinalVoxInds,:);
                end
            end
            
            clear sorted* 
        end
    end
    
    % Now plot summaries for this model (official min number of vox = 75)
    
    
end
end

%% Make null tval distributions for each ROI and featType

% This probably goes after next part but before plot!


nPerms = size(nullRSqArrayFinalVoxEachRoiVarTypeAndSubj{1,1},2);
nullTvalsEachRoiAndVarType = ...
    nan(numel(roiNames), numel(varTypeLabels), nPerms);
tvalsEachRoiAndVarType = nan(numel(roiNames), numel(varTypeLabels));
pvalsEachRoiAndVarType = nan(numel(roiNames), numel(varTypeLabels));

nVarCombos = 3; % 3 var types, 3 pairwise combos.
nullTvalsDiffsEachRoiAndVarCombo = ...
    nan(numel(roiNames), nVarCombos, nPerms);
tvalsDiffsEachRoiAndVarCombo = nan(numel(roiNames), nVarCombos);
pvalsDiffsEachRoiAndVarCombo = nan(numel(roiNames), nVarCombos);


for iRoi = 1:numel(roiNamesToUse)
    
    for iVarType = 1:numel(varTypeLabels)
        curNullMeansAcrossVoxBySubjAndPerm = ...
            nan(numel(subjNumsStr), nPerms); % this is only for the current roi and var type 
        for iSubj = 1:numel(subjNumsStr)
            % Have to mean across vox to get one val per roi, subj, and
            % perm (and then tstats will be done across these vals across
            % subjs)
            if roiIsUsableBySubjAndRoi(iSubj,iRoi)
                
                % automated error if taking mean across pre-vox-cutoff
                % voxels
                if size(nullRSqArrayFinalVoxEachRoiVarTypeAndSubj{iRoi,iVarType,iSubj},1) > maxVoxToUse
                   error('Really bad error! You''re supposed to be mean-ing for each perm AFTER applying the vox cutoff');
                end
                curNullMeansAcrossVoxBySubjAndPerm(iSubj,:) = ...
                    mean(nullRSqArrayFinalVoxEachRoiVarTypeAndSubj{iRoi,iVarType, iSubj}); % inside each cell is just vox; in contrast to regular analysis, here the var type is a cell dim and not inside the cell
            else
                % leave NaNs here to count in lines below if there aren't
                % enough vox this ROI/subj.
            end
        end
        
        curSubjWEnoughVox = unique(sum(~isnan(curNullMeansAcrossVoxBySubjAndPerm),1)); % this is the number of subjects with enough vox (should be 7 or 8)
        if numel(curSubjWEnoughVox) > 1 % This is to say how many of the 8 subjects have enough voxels for this ROI. Should be 8 for all except RSC, which should have 7.
            error('Oh no! Different numbers of usable subjs for diff vox. That''s an error!'); % (would error later anyway, but this way won't forget why there'd be an error here.)
        end
        
        % get tstat across subjs
        % (remember I have to only take subjects with >= 75 vox though!)
        % getting one for each random permutation, to get a whole
        % distribution.
        nullTvalsEachRoiAndVarType(iRoi,iVarType,:) = ...
            (nanmean(curNullMeansAcrossVoxBySubjAndPerm,1) - 0) ./ ... % comparing against 0 for each perm
            (nanstd(curNullMeansAcrossVoxBySubjAndPerm,1) ./ ...
            sqrt(curSubjWEnoughVox));
        
        % Now get real (non-null) tvals for this ROI and var type (bc
        % that's what we're doing each "t-test" within)
        curRSqs = squeeze(rSqToPlotEachRoiVarTypeAndSubj(iRoi, iVarType, :)); % this line needs to be double-checked to get correct array name.
        tvalsEachRoiAndVarType(iRoi,iVarType) = ...
            nanmean(curRSqs) ./ (nanstd(curRSqs) ./ sqrt(curSubjWEnoughVox));
        
        % Get number of perms more extreme and also pvals.
        pvalsEachRoiAndVarType(iRoi,iVarType) = getTwoTailedPermPval( ...
            tvalsEachRoiAndVarType(iRoi,iVarType), ...
            nullTvalsEachRoiAndVarType(iRoi,iVarType,:));
        
    end
    
    
    % Now get tvals for differences comparing variance types
    % ------------------------------------------------------------
    varComboStrings = cell(3,1);
    for iVarCombo = 1:nVarCombos
        
        varComboMatrix = [ ... % All possible comparisons of the 3 var types (BU, shared, AU)
            1 2;
            1 3;
            2 3];
        curVarCombo = varComboMatrix(iVarCombo,:);
        varComboStrings{iVarCombo} = ...
            [varTypeLabels{curVarCombo(1)} ' vs ' varTypeLabels{curVarCombo(2)}];
        
        curSubtractedNullValsBySubjAndPerm  = nan(numel(subjNumsStr), nPerms);
        for iSubj = 1:numel(subjNumsStr)
            % Tvals from subtracted vals (operations done across subjs)
            curSubtractedNullValsBySubjAndPerm(iSubj,:) = ... % Taking mean (across vox in this case) and then subtracting is the same as subtracting and then taking the mean
                mean(nullRSqArrayFinalVoxEachRoiVarTypeAndSubj{iRoi,curVarCombo(1), iSubj},1) - ...
                mean(nullRSqArrayFinalVoxEachRoiVarTypeAndSubj{iRoi,curVarCombo(2), iSubj},1); % End up w one value
        end
        
        nullTvalsDiffsEachRoiAndVarCombo(iRoi,iVarCombo,:) = ... % Get a paired-difference tval for each perm in null distribution
            nanmean(curSubtractedNullValsBySubjAndPerm,1) ./ ...
            (nanstd(curSubtractedNullValsBySubjAndPerm,1) ./ sqrt(curSubjWEnoughVox));
        
        % Now get real r-squared differences and turn into tvals.
        curRSqsDiffs = ...
            squeeze(rSqToPlotEachRoiVarTypeAndSubj(iRoi, curVarCombo(1), :)) - ...
            squeeze(rSqToPlotEachRoiVarTypeAndSubj(iRoi, curVarCombo(2), :)); % this line needs to be double-checked to get correct array name.
        
        tvalsDiffsEachRoiAndVarCombo(iRoi, iVarCombo) = ...
            nanmean(curRSqsDiffs,1) ./ ...
            (nanstd(curRSqsDiffs,1) ./ sqrt(curSubjWEnoughVox));
        
        % Now take real tvals and null-dist of tvals and compare to get
        % pvals.
        pvalsDiffsEachRoiAndVarCombo(iRoi, iVarCombo) = getTwoTailedPermPval( ...
            tvalsDiffsEachRoiAndVarCombo(iRoi,iVarCombo), ...
            nullTvalsDiffsEachRoiAndVarCombo(iRoi,iVarCombo,:));
    end
end
%}

% Now make/save pvals and stats
% -----------------------------------------------------

individVarTypePvals = array2table(pvalsEachRoiAndVarType, ...
    'VariableNames', varTypeLabels, ...
    'RowNames', roiNamesToUse);
disp(individVarTypePvals);

varComparisonPvals = array2table(pvalsDiffsEachRoiAndVarCombo, ...
    'VariableNames', genvarname(varComboStrings), ...
    'RowNames', roiNamesToUse);
disp(varComparisonPvals);


%% Now plot summaries and also plot stats on top!


if doAvgSignalAnalysis && plotThresholdVoiDetails
    
    %{
    for iSubj = subjsToDoNow
        
        curSubj = subjIds{iSubj};
        
        
        % figure; 
        % plotCtr = 0;
        
        for iRoi = roisToDoNow
            curRoi = roiNamesToUse{iRoi};
            
            % plotCtr = plotCtr + 1;
            % subplot(2, 3, plotCtr); hold on; % each ROI gets its own pane.
            
            
            for iVarType = 1:numel(varTypeLabels)
                
                % featPrefix = ['_' featLabels{iFeatToPlot} '_'];

                
                % load intermediate file! 
                % ------------------------------
                
                % Summaries of prediction performances (r's)
                % Specific *just* to feature type, but not taking the trouble to only read in once.
                %{
                curRsIntermedFName = fullfile(intermedDir, ...
                    ['acrossSubjRoiSummaries' featPrefix precText ...
                    versionText hemiText intText '_' trWindowText '_' ...
                    peakTrText cvFoldText sceneText '.mat']);
                load(curRsIntermedFName, ...
                    'rsBySubjAndRoi', 'meanRsBySubjAndRoi', 'medRsBySubjAndRoi', ...
                    'semRsBySubjAndRoi', 'subjsToDoNow', 'roiNamesToUse'); % plotting meanRsByVox, also using numCvSplitsOuter to scale the noise ceiling estimates to match a corr w the same number of timepoints as these predictions
                
                
                % Summaries of localizer activations
                % Specific to feature type, subj, ROI
                curLocTvalsIntermedFName = fullfile(intermedDir, [ ...
                    'rSqsSortedByVoxLocActivation' curSubj '_' curRoi precText ...
                    versionText hemiText intText '_' trWindowText '_' ...
                    sceneText '.mat']);
                load(curLocTvalsIntermedFName, ...
                    'curSubj','curSubjNum','iSubj','curArrayOfAvgsByVarType', ...
                    'sortedEvcTvals', 'highestEvcTvalInds', ...
                    'sortedSceneTvals', 'highestSceneTvalInds', ...
                    'varTypeLabels','iRoi','curRoi'); % plotting meanRsByVox, also using numCvSplitsOuter to scale the noise ceiling estimates to match a corr w the same number of timepoints as these predictions
                %}
                
                
            end
            % legend(plotHandle, featLabels);
            
        end
    end
    %}
    if featTypeA == 1
        varTypeBarColors = [ ...
            .7 .2 .6; ...
            .7 .7 .65; ...
            .7 .7 .2];
        varTypeScatterColors = [ ...
            .5 .1 .4; ...
            .5 .5 .4; ...
            .4 .4 .1];
    elseif featTypeA == 2
        varTypeBarColors = [ ...
            .7 .2 .6; ...
            .7 .7 .65; ...
            .8 .5 .2];
        varTypeScatterColors = [ ...
            .5 .1 .4; ...
            .5 .5 .4; ...
            .6 .4 .1];
    end
    voiFontSize = 12;
    axLabelFontSize = 12;
    
    % if ~doFittingByTask
        if ~exist('rSqToPlotEachRoiVarTypeAndSubj','var') 
            error('Woops! You need to run the previous section first.')
        end
        plotThingsOnTopOfGroupedBarPlot(rSqToPlotEachRoiVarTypeAndSubj, 3, 1, 2, ...
            pvalsEachRoiAndVarType, varTypeBarColors, varTypeScatterColors);
        title({'Var. explained', 'up to 200 best localizer voxels','min 75'});
        ylabel('Propor variance explained, r-squared');
        % xlabel('ROI'); 
        xticklabels(roiNamesToUse); % , 'fontsize', voiFontSize);
        curAx = gca; curAx.XAxis.FontSize = voiFontSize; 
        curAx.XAxisLocation = 'origin'; curAx.XAxis.TickLength = [0 0];
        % curAx.XAxis.FontWeight = 'bold'; curAx.XLabel.FontWeight = 'normal';
        curAx.XLabel.FontSize = axLabelFontSize; 
        curAx.YAxis.FontSize = axLabelFontSize; % (this does label and ticks this size)
    % else
        % error('oh no! variance partitioning isn''t set up for task comparisons yet.');
    % end
    
    legend(varTypeLabels, 'FontSize', axLabelFontSize);
    
    
    
    %{
    groupedBarWErrorBarsSubjDots(rSqToPlotEachRoiVarTypeAndSubj, ...
        roiNamesToUse, newRoiOrders, ...
        grays, ...
        curXLabel, curYLabel, curLegendLabels, curTitle);
    %}
    
end 




%% Plot lambdas, also plot against performances
% (This part only works for the models fitted in Matlab.)
%{
if ~fromPython
    for iSubj = 1:numel(subjIds)
        curSubj = subjIds{iSubj};
        fprintf('Starting %s \n',curSubj);
        % figure;
        
        for iRoi = 1:numel(roiNamesToUse)
            curRoi = roiNamesToUse{iRoi};
            curFName = fullfile( ...
                featDir, ['fittedModels' featPrefix curSubj '_' curRoi precText versionText hemiText '.mat']);
            
            load(curFName, 'lambdasBySplit');
            
            figure; semilogy(lambdasBySplit);
            keyboard;
        end
    end
end
%}

%% Plots

function plotNullDistAndStats(permArray,realValuesByVox, ...
    pairedColorMap,colorIndsOldOrders,iRoi,curRoi,curSubj,varTypeLabel)

% Now do null distribution of vox means figure for this ROI, print out stats.
figure; hold on;

curNullVoxMeanArray = mean(permArray,1); % Mean across rows, which are vox; leaves 10,000 null-dist means.
curMean = mean(realValuesByVox); % also taking mean across voxels, leaving one value.
permsBelowMean = sum(curNullVoxMeanArray < curMean);
lessExtremePerms = sum(abs(curNullVoxMeanArray) < abs(curMean));
nPerms = size(permArray,2);
curPval = 1 - ((lessExtremePerms-1)/ nPerms);
curPrctile975 = prctile(curNullVoxMeanArray,97.5);
% prctile975BySubjAndRoi(iSubj,iRoi) = curPrctile975;

histogram(curNullVoxMeanArray,'FaceColor',[0.3 0.3 0.3], 'EdgeColor', [.4 .4 .4]);
line([curMean curMean], ylim, 'Color', pairedColorMap(colorIndsOldOrders(iRoi),:), ...
    'LineWidth', 3, 'LineStyle', ':');
line([curPrctile975 curPrctile975], ylim, 'Color', [0.3 0.3 0.3], ...
    'LineStyle',':');
curYLim = ylim;
text(curMean - (curMean/50), curYLim(2)-(curYLim(2)/25), ...
    sprintf('%d/%d permutations\n below actual mean;\nTwo-tailed p < %0.4f', ...
    permsBelowMean, nPerms, curPval), 'FontWeight','bold', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment','top');
title(sprintf('Null distribution of vox means,\n%s, %s, %s', ...
    varTypeLabel, curSubj,curRoi));
xlabel('Var (r^2), mean across vox');
end

function groupedBarWErrorBarsSubjDots(arrayOfValsInOrder, ...
    roiNamesToUse, newRoiOrders, ...
    barColors, ...
    curXLabel, curYLabel, curLegendLabels, curTitle)

% For each ROI, separately get within-subjs sem, normalizing by mean of var
% explained across variance types
semValsToPlot = nan(size(arrayOfValsInOrder,1), size(arrayOfValsInOrder,2));
for iRoi = 1:numel(roiNamesToUse)
    semValsToPlot(iRoi,:) = ... % gives val for each variance type
        semWithin(squeeze(arrayOfValsInOrder(iRoi,:,:))'); % semwithin wants subjects x conditions, so transpose after squeeze for each ROI's data.
end

figure;

box off; hold on;

meanValsToPlot = mean(arrayOfValsInOrder,3); % subj should be third dim!
b = bar(meanValsToPlot(newRoiOrders,:)); % ROIs should be first dim, then var type; will error if dims are wrong
xticklabels(roiNamesToUse(newRoiOrders));


ngroups = size(meanValsToPlot,1);
nbars = size(meanValsToPlot, 2);
groupwidth = min(.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange

for iBar = 1:nbars % Bar is probably variance type (A-unique, B-unique, shared) but could just as easily be any condition label if plotting something else
    % Calculate center of each bar
    % (Will be able to tell if sem is wrong, bc order should be wrong for
    % both)
    x = (1:ngroups) - groupwidth/2 + (2*iBar-1) * groupwidth / (2*nbars);
    errorbar(x, meanValsToPlot(newRoiOrders,iBar), ...
        semValsToPlot(newRoiOrders,iBar), 'Color',[0.2 0.2 0.2], ...
        'CapSize',3,'linestyle', 'none'); % 'Color',pairedColorMap(i,:),
    
    for iSubj = 1:size(arrayOfValsInOrder,3) % 3rd dim is subj
        scatter(x, squeeze(arrayOfValsInOrder(:,iBar,iSubj)), 16, [0.1 0.1 0.1], ...
            '.')
    end
    
    b(iBar).FaceColor = barColors(iBar,:);
end


% legend(roiNamesToUse(newRoiOrders));
% ylabel('Mean of corrs w held-out data');
xlabel(curXLabel);
% xticklabels(roiNames(newRoiOrders));
set(gca,'xticklabel',{[]})
ylabel(curYLabel); % 'Propor variance explained, r-squared');
legend(b,curLegendLabels);
title(curTitle);
end

%% Functions for sorting VOI voxels by independent-localizer activation 

function [sortedEvcTvals, highestEvcTvalInds, ...
    sortedSceneTvals, highestSceneTvalInds] = ...
    getSortedTvalsFromArrays(curEvcLocTvals, curSceneLocTvals)

% Inputs are evc loc tvals that are already for the current subject and
% ROI. 

[sortedEvcTvals, highestEvcTvalInds] = ...
                sort(curEvcLocTvals, 'descend');
[sortedSceneTvals, highestSceneTvalInds] = ...
                sort(curSceneLocTvals, 'descend');
            
end

function [curArrayOfAvgsByVarType, plotHandle] = plotVoxSelectionStackedPlotWithinLoop(rSqsBySubjAndRoiAndVarType, ...
    sortedCorrespondingTvals, highestCorrespondingTvalInds, ... % these two inputs should be controlled outside loop (either evc or scene tvals depending on ROI).
    varTypeLabels, iRoi, iSubj, curRoi, curSubjNum, plotCtr)

locLabelsByRoi = {'evc', 'evc', 'evc', 'scene', 'scene', 'scene'};
% Figure showing effects of vox selection
% -------------------------------------------------------

subplot(2,3,plotCtr); hold on;

nVarTypes = size(rSqsBySubjAndRoiAndVarType{iSubj,iRoi},2); % column of each cell will be var type; probably 3 of them.
curArrayOfAvgsByVarType = nan(numel(sortedCorrespondingTvals),nVarTypes); % Size of number of voxels bc starting by including one, ending by including all

for iNHighest = 1:numel(sortedCorrespondingTvals)
    voxToInclude = highestCorrespondingTvalInds(1:iNHighest);
    for iVarType = 1:nVarTypes % Keeps var types in whatever order they were in
        curArrayOfAvgsByVarType(iNHighest,iVarType) = ...
            mean(rSqsBySubjAndRoiAndVarType{iSubj,iRoi}(voxToInclude, iVarType));
    end
end
plotHandle = area(curArrayOfAvgsByVarType);
xlabel(['Top N ' locLabelsByRoi{iRoi} ' loc. voxels included in avg. perf.']);


ylabel('Variance explained (r-squared)');
curYlim = ylim;
line([75 75], ylim, 'Color', [0.8 0.8 0.8], 'LineStyle', ':');
line([200 200], ylim, 'Color', [0.3 0.3 0.3]); % max 200 vox used, even though more in arrays.
xlim([-0.5 200.5]); % 4.5]);
legend(varTypeLabels); % Making legend the last thing drawn keeps it from including extra entries
title(sprintf('%s, %s', curSubjNum, curRoi));

end
