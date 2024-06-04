% 2024.06.02 AS updating visualization script to be cleaner/clearer to output figures

%%
clear;
% close all;

addpath(genpath('helpers')); % Add path to helper function
subjIds = {'CI','CO','CK','CR','BX','CS','CG','CH'}; 
subjNumsStr = {'S01','S02','S03','S04','S05','S06','S07','S08'}; 
subjsToDoNow = 1:8;
roiNames = {'V1','V2','V3', ...
    'OPA', 'PPA', 'RSC'};
roiToStart = 1;
roiToEnd = 6;

doVoxSelection = 1; % This selects the top N voxels based on independent localizer data (how data is analyzed for main analysis)
doMainAnalysis = 1; % Do main analysis and visualization (Fig3A; not just saving intermediate outputs for later)
saveIntermedFiles = true; % To save time, the outputs of the analyses done in this script are saved into intermediate files! Decide here if you want to create these intermediates.
useSavedIntermedFiles = true; % To save time, you can plot figures using files saved from a previous run of this script! You can set this to "false" if you don't want to use the previously generated files.
doStats = true; % This does permutation tests, which can take awhile, so you can set to "false" if you just want to see the participants' data without stats.
numCvSplitsOuter = 9; % 9; Can be 18 (2 runs grouped together) or 9 (4 runs grouped together); there are 36 runs total.

analysisDir = pwd;
filesepinds = find(analysisDir==filesep);
mainDir = analysisDir(1:filesepinds(end)-1);
fmriDataDir = fullfile(mainDir,'DataFmri');
intermedDir = fullfile(analysisDir,'ridgeFigs','intermedArraysAcrossModels');

trWindowText = '2.4-4.8sec'; % Suffix for the window of TRs used for the main analysis

if numCvSplitsOuter == 9
    cvFoldText = 'cv9';
elseif numCvSplitsOuter == 18
    cvFoldText = 'cv18';
end

hemiText = '_collapsedHemispheres';
roiNamesToUse = roiNames;

% Putting together all suffixes into one string of options to add to file names
analysisOptionsText = [hemiText '_' trWindowText '_' cvFoldText]; 


% ------------------------------------
% Feature-set info/options
% ------------------------------------

featLabels = {'spat3d','nonSpat3d','gabor'}; % ,'pix'}; % This will keep pix feats from being plotted.
featTypesToMakeIntermedFilesFor = 1:3; % Can do all of these here or


%% First portion of script loops through each feature type
for featType = featTypesToMakeIntermedFilesFor
    
    featPrefix = ['_' featLabels{featType} '_'];
    featDir = fullfile('ridgeOutputs', ...
        [featLabels{featType} 'Feats' ...
        '_' trWindowText '_' cvFoldText]);
    disp(['Starting ' featLabels{featType} ' features...']);
    
    % ------------------------------------
    % Variable initialization for this feature type
    % ------------------------------------
    
    % These arrays include all voxels, and they'll be selected based on
    % independent localizer data afterwards.
    rsBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse)); % also includes voxels inside each cell
    if doStats
        nullRsBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse));
    end
    
    
    %% Make/fill summaries; plot hists for each ROI/subj
    % this is to grab the individual outputs of the model fitting and put
    % it together in arrays so it's easier to plot
    % Also grabs the localizer activations that will later be used for
    % voxel selection.
    
    voxDistsFig = figure; hold on; % For all distributions of fitting success across voxels
    for iSubj = subjsToDoNow
        curSubj = subjIds{iSubj};
        curSubjNum =  subjNumsStr{iSubj};
        fprintf('   Starting %s \n',curSubj);
        
        for iRoi = roiToStart:roiToEnd 
            curRoi = roiNamesToUse{iRoi};
            
            if doMainAnalysis % The input text to the file name (controlled above) will determine which gets done
                
                curInFName = fullfile( ...
                    featDir, ['fittedModels' featPrefix curSubj '_' curRoi ...
                    analysisOptionsText '.mat']);
                
                load(curInFName, 'meanRsByVox', 'numCvSplitsOuter'); % plotting meanRsByVox, also optionally using numCvSplitsOuter to scale the noise ceiling estimates to match a corr w the same number of timepoints as these predictions
                
                if doStats
                    load(curInFName,'meanRsByVox_perms','numCvSplitsOuter');
                    nullRsBySubjAndRoi{iSubj,iRoi} = meanRsByVox_perms;
                    if iSubj==1 && iRoi==1
                        nPerms = size(meanRsByVox_perms,2);
                    end
                end
                
                % Each voxel's r from model fitting
                rsBySubjAndRoi{iSubj,iRoi} = meanRsByVox; % This variable is called "mean..." bc it's the mean across cv folds! Will take mean across vox after vox selection.
                
                clear meanRsByVox*
            end
            
        end; clear curRoi;
    end; clear curSubj curSubjNum;
    
    % ------------------------------------------------------------------------
    % Save intermediate summaries into files (rsBySubjAndRoi)
    % so that I can plot and do stats across models in later sections.
    % ------------------------------------------------------------------------
    
    if saveIntermedFiles
        
        if doMainAnalysis 
            
            curIntermedFName = ...
                fullfile(intermedDir, ...
                ['acrossSubjRoiSummaries' featPrefix analysisOptionsText '.mat']);
            
            if doStats
                save(curIntermedFName, 'nullRsBySubjAndRoi', ...
                    'rsBySubjAndRoi', ... 
                    'subjsToDoNow', 'roiNamesToUse', ...
                    'featType','featLabels', '-v7.3');
            else
                save(curIntermedFName, ...
                    'rsBySubjAndRoi', ...
                    'subjsToDoNow', 'roiNamesToUse', ...
                    'featType','featLabels');
            end
            clear curIntermedFName
        end
    else
    end
    
    
    % ------------------------------------------------------------------------
    % Vox included in order of loc activity;
    % This section only plots individual subjects/models at a time (saving info
    % in intermediate files to be able to plot all of them later)
    % ------------------------------------------------------------------------
    
    if doVoxSelection
        
        % Selecting voxels based on independent localizer data
        % (Reading in t-values from the contrasts in the next 2 lines)
        sceneLocTvalsBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse)); % scenes > objects
        evcLocTvalsBySubjAndRoi = cell(numel(subjIds),numel(roiNamesToUse)); % scrambled > baseline
        
        for iSubj = 1:numel(subjIds)
            voxSelectionPlot = figure;
            plotCtr = 0;
            curSubjNum = subjNumsStr{iSubj};
            curSubj = subjIds{iSubj};
            
            for iRoi = 1:numel(roiNamesToUse)
                curRoi = roiNamesToUse{iRoi};
                plotCtr = plotCtr + 1;
                
                % File is labeled "main task", but only reading in the
                % localizer info!
                curLocTvalFName = fullfile(fmriDataDir, ...
                    ['MainTaskSignalByTrial_' trWindowText ...
                    '_' curSubjNum, '.mat']);
                load(curLocTvalFName, 'localizerInfo');
                
                % Info about localizer t-vals (for voxel selection)
                % There's info saved for both contrasts for each ROI, but
                % we'll only ultimately end up using scenes > objects for
                % the scene regions and scrambled > baseline for EVC
                % regions (see below)
                sceneLocTvalsBySubjAndRoi{iSubj,iRoi} = localizerInfo(iRoi).sceneTvals; % scenes > objects
                evcLocTvalsBySubjAndRoi{iSubj,iRoi} = localizerInfo(iRoi).evcTvals; % scrambled > baseline
                
                
                [sortedEvcTvals, highestEvcTvalInds, ...
                    sortedSceneTvals, highestSceneTvalInds] = ...
                    getSortedTvalsFromArrays( ...
                    evcLocTvalsBySubjAndRoi{iSubj,iRoi}, ...
                    sceneLocTvalsBySubjAndRoi{iSubj,iRoi});
                
                
                % Figures showing how/whether effects vary depending on number
                % of voxels included (just one feat type here)
                % -------------------------------------------------------
                
                set(0,'CurrentFigure',voxSelectionPlot);
                if iRoi < 4 % Use evc tvals as inputs for EVC ROIs
                    [curArrayOfAvgs, ~] = plotVoxSelectionPlotWithinLoop(rsBySubjAndRoi, ...
                        sortedEvcTvals, highestEvcTvalInds, ...
                        iRoi, iSubj, curRoi, curSubjNum, plotCtr);
                else % Use scene tvals as inputs for scene-selective ROIs
                    [curArrayOfAvgs, ~] = plotVoxSelectionPlotWithinLoop(rsBySubjAndRoi, ...
                        sortedSceneTvals, highestSceneTvalInds, ...
                        iRoi, iSubj, curRoi, curSubjNum, plotCtr);
                end
                
                if saveIntermedFiles
                    
                    % To save the files with voxels sorted by their respective localizer activations
                    % (Used for main analyses where choosing the top
                    % 200 voxels!)
                    curIntermedFName = ...
                        fullfile(intermedDir, ...
                        ['rsSortedByVoxLocActivation' featPrefix curSubj '_' curRoi ...
                        analysisOptionsText '.mat']);
                    
                    save(curIntermedFName, ...
                        'curSubj','curSubjNum','iSubj','curArrayOfAvgs', ...
                        'featType','featLabels','iRoi','curRoi', ...
                        'sortedEvcTvals', 'highestEvcTvalInds', ...
                        'sortedSceneTvals', 'highestSceneTvalInds');
                    clear curIntermedFName
                end
                
                clear sorted*
            end; clear curRoi;
        end; clear curRoi;
        
    end
end

%% If saved for all feature types, plot on top of each other!

if useSavedIntermedFiles && doMainAnalysis && doVoxSelection
    
    % This is for the plot that only uses the top 200 voxels (so the main plot).
    rToPlotEachRoiModelAndSubj = ...
        nan(numel(roiNamesToUse), numel(featLabels), numel(subjIds));
    nullRArrayFinalVoxEachRoiModelAndSubj = ...
        cell(numel(roiNamesToUse), numel(featLabels), numel(subjIds));
    
    % This is a corresponding array to say how many voxels used
    nVoxBySubjAndRoi = nan(numel(subjIds), numel(roiNamesToUse));
    roiIsUsableBySubjAndRoi = nan(size(nVoxBySubjAndRoi));
    
    for iSubj = subjsToDoNow
        
        curSubj = subjIds{iSubj};
        
        figure;
        plotCtr = 0;
        
        for iRoi = roiToStart:roiToEnd
            curRoi = roiNamesToUse{iRoi};
            
            plotCtr = plotCtr + 1;
            subplot(2, 3, plotCtr); hold on; % each ROI gets its own pane.
            
            
            for iFeatToPlot = 1:numel(featLabels)
                
                featPrefix = ['_' featLabels{iFeatToPlot} '_'];
                
                % load intermediate file!
                % ------------------------------
                
                curRsIntermedFName = fullfile(intermedDir, ...
                    ['acrossSubjRoiSummaries' featPrefix ...
                    analysisOptionsText '.mat']); % Note: hard-coding for debugging! AS change later!% % sceneText '.mat']);
                
                load(curRsIntermedFName, ...
                    'rsBySubjAndRoi', ... 'meanRsBySubjAndRoi', ...'medRsBySubjAndRoi', 'semRsBySubjAndRoi', ...
                    'subjsToDoNow', 'roiNamesToUse'); % plotting meanRsByVox, also using numCvSplitsOuter to scale the noise ceiling estimates to match a corr w the same number of timepoints as these predictions
                if doStats
                    load(curRsIntermedFName, 'nullRsBySubjAndRoi');
                end
                
                % Summaries of localizer activations
                % Specific to feature type, subj, ROI
                % (Sorted by localizer activation to be able to plot
                % the impact of including different numbers of voxels)
                
                curLocTvalsIntermedFName = fullfile(intermedDir, [ ...
                    'rsSortedByVoxLocActivation' featPrefix curSubj '_' curRoi ...
                    analysisOptionsText '.mat']); % Note: hard-coding for debugging! AS change later!% % sceneText '.mat']);
                
                load(curLocTvalsIntermedFName, ...
                    'curSubj','curSubjNum','iSubj','curArrayOfAvgs', ...
                    'sortedEvcTvals', 'highestEvcTvalInds', ...
                    'sortedSceneTvals', 'highestSceneTvalInds', ...
                    'featType','iRoi','curRoi'); % plotting meanRsByVox, also using numCvSplitsOuter to scale the noise ceiling estimates to match a corr w the same number of timepoints as these predictions
                
                % Basic info about vox selection
                % -------------------------------------------------------
                % Save info about how many vox used here, as well as the
                % r-value corresponding to that number of vox.
                % (Redundant to do some of this for each model, but the
                % info's only read in from files saved for each model, so
                % just doing it redundantly)
                maxVoxToUse = 200; % Max vox to average over
                voxCutoff = 75; % Min vox to be included in analysis.
                curNVoxUsed = min(maxVoxToUse, size(curArrayOfAvgs,1)); % Either use number of vox (if less than 200) or use 200 (cutoff)
                nVoxBySubjAndRoi(iSubj,iRoi) = curNVoxUsed;
                roiIsUsableBySubjAndRoi(iSubj,iRoi) = curNVoxUsed >= voxCutoff;
                if roiIsUsableBySubjAndRoi(iSubj,iRoi) % Fill in if usable; otherwise leave NaN
                    rToPlotEachRoiModelAndSubj(iRoi,iFeatToPlot,iSubj) = ...
                        curArrayOfAvgs(curNVoxUsed); % curNVoxUsed is the number of voxels to average over, and curArrayOfAvgs is ordered by how many vox to avg over.
                end
                
                % Figure showing effects of vox selection 
                % -------------------------------------------------------
                
                if iRoi < 4
                    [~, plotHandle(iFeatToPlot)] = plotVoxSelectionPlotWithinLoop(rsBySubjAndRoi, ...
                        sortedEvcTvals, highestEvcTvalInds, ... % Use evc tvals as inputs for EVC ROIs
                        iRoi, iSubj, curRoi, curSubjNum, plotCtr);
                    curFinalVoxInds = highestEvcTvalInds(1:curNVoxUsed);
                else
                    [~, plotHandle(iFeatToPlot)] = plotVoxSelectionPlotWithinLoop(rsBySubjAndRoi, ...
                        sortedSceneTvals, highestSceneTvalInds, ... % Use scene tvals as inputs for scene ROIs
                        iRoi, iSubj, curRoi, curSubjNum, plotCtr);
                    curFinalVoxInds = highestSceneTvalInds(1:curNVoxUsed);
                end
                
                
                % NOW USE THIS TO GET BY-SUBJ NULL VALS EACH FEAT SET
                % AND ROI
                % -----------------------------------------------------
                % Can use highestScene(/Evc)TvalInds to get top N
                % voxels here, to save correct vox for null distribution.
                % (also check method indexing into rsBy... vs.
                % method using curArrayOfAvgs...
                
                if doStats
                    nullRArrayFinalVoxEachRoiModelAndSubj{iRoi,iFeatToPlot,iSubj} = ...
                        nullRsBySubjAndRoi{iSubj,iRoi}(curFinalVoxInds,:);
                end
                
            end
            legend(plotHandle, featLabels);
            
        end; clear curRoi;
    end; clear curSubj curSubjNum;
    
end

%% Make null tval distributions for each ROI and featType

if useSavedIntermedFiles && doStats
    nullTvalsEachRoiAndModel = ...
        nan(numel(roiNames), numel(featLabels), nPerms);
    tvalsEachRoiAndModel = nan(numel(roiNames), numel(featLabels));
    pvalsEachRoiAndModel = nan(numel(roiNames), numel(featLabels));
    
    nModelCombos = 3; % 3 model types, 3 pairwise combos.
    nullTvalsDiffsEachRoiAndModelCombo = ...
        nan(numel(roiNames), nModelCombos, nPerms);
    tvalsDiffsEachRoiAndModelCombo = nan(numel(roiNames), nModelCombos);
    pvalsDiffsEachRoiAndModelCombo = nan(numel(roiNames), nModelCombos);
    
    
    for iRoi = roiToStart:roiToEnd 
        
        for iFeatType = 1:numel(featLabels)
            curNullMeansAcrossVoxBySubjAndPerm = ...
                nan(numel(subjNumsStr), nPerms);
            for iSubj = subjsToDoNow 
                
                % Have to mean across vox (after vox selection) to get one 
                % val per roi, subj, and perm (and then tstats will be done 
                % across these vals across subjs)
                if roiIsUsableBySubjAndRoi(iSubj,iRoi)
                    
                    % Double-check the voxel cutoff has already been
                    % applied
                    if size(nullRArrayFinalVoxEachRoiModelAndSubj{iRoi,iFeatType,iSubj},1) > maxVoxToUse
                        error('Really bad error! You''re supposed to be mean-ing for each perm AFTER applying the vox cutoff');
                    end
                    
                    % Also for null arrays, take mean across vox (1st
                    % dim)
                    curNullMeansAcrossVoxBySubjAndPerm(iSubj,:) = ...
                        mean(nullRArrayFinalVoxEachRoiModelAndSubj{iRoi,iFeatType,iSubj},1); % first dim is vox; 2nd is perm
                else
                    % If an ROI isn't usable, this empty "else" will leave 
                    % NaNs here to count in lines below if there aren't
                    % enough vox this ROI/subj.
                end
            end
            
            % This portion is to check how many of the 8 subjects have enough voxels for this ROI. Should be 8 for all except RSC, which should have 7.
            curSubjWEnoughVox = unique(sum(~isnan(curNullMeansAcrossVoxBySubjAndPerm),1));
            if numel(curSubjWEnoughVox) > 1 
                error('Oh no! Different numbers of usable subjs for diff vox. That''s an error!'); % (would error later anyway, but this way won't forget why there'd be an error here.)
            end
            
            % Get tstat across subjs to use as test statistic (doing
            % non-parametric tests on the tstat though!)
            % (remember to only take subjects with >= 75 vox)
            nullTvalsEachRoiAndModel(iRoi,iFeatType,:) = ...
                (nanmean(curNullMeansAcrossVoxBySubjAndPerm,1) - 0) ./ ... % comparing against 0 for each perm
                (nanstd(curNullMeansAcrossVoxBySubjAndPerm,1) ./ ...
                sqrt(curSubjWEnoughVox));
            
            % Now get real tvals
            curRVals = squeeze(rToPlotEachRoiModelAndSubj(iRoi, iFeatType, :));
            tvalsEachRoiAndModel(iRoi,iFeatType) = ...
                nanmean(curRVals) ./ (nanstd(curRVals) ./ sqrt(curSubjWEnoughVox));
            
            % Get number of perms more extreme and also pvals.
            pvalsEachRoiAndModel(iRoi,iFeatType) = getTwoTailedPermPval( ...
                tvalsEachRoiAndModel(iRoi,iFeatType), ...
                nullTvalsEachRoiAndModel(iRoi,iFeatType,:));
        end
        
        
        % For same ROI, now get tvals for differences comparing variance types
        % ------------------------------------------------------------
        modelComboStrings = cell(3,1);
        for iModelCombo = 1:nModelCombos
            
            modelComboMatrix = [ ... % All possible comparisons of the 3 var types (BU, shared, AU)
                1 2;
                1 3;
                2 3];
            curModelCombo = modelComboMatrix(iModelCombo,:);
            modelComboStrings{iModelCombo} = ...
                [featLabels{curModelCombo(1)} ' vs ' featLabels{curModelCombo(2)}];
            
            curSubtractedNullValsBySubjAndPerm  = nan(numel(subjNumsStr), nPerms);
            for iSubj = subjsToDoNow 
                % Tvals from subtracted vals (operations done across subjs)
                curSubtractedNullValsBySubjAndPerm(iSubj,:) = ... % Taking mean (across vox in this case) and then subtracting is the same as subtracting and then taking the mean
                    mean(nullRArrayFinalVoxEachRoiModelAndSubj{iRoi,curModelCombo(1), iSubj},1) - ...
                    mean(nullRArrayFinalVoxEachRoiModelAndSubj{iRoi,curModelCombo(2), iSubj},1); % End up w one value
            end
            
            nullTvalsDiffsEachRoiAndModelCombo(iRoi,iModelCombo,:) = ... % Get a paired-difference tval for each perm in null distribution
                nanmean(curSubtractedNullValsBySubjAndPerm,1) ./ ...
                (nanstd(curSubtractedNullValsBySubjAndPerm,1) ./ sqrt(curSubjWEnoughVox));
            
            curRsDiffs = ...
                squeeze(rToPlotEachRoiModelAndSubj(iRoi, curModelCombo(1), :)) - ...
                squeeze(rToPlotEachRoiModelAndSubj(iRoi, curModelCombo(2), :)); % this line needs to be double-checked to get correct array name.
            
            tvalsDiffsEachRoiAndModelCombo(iRoi, iModelCombo) = ...
                nanmean(curRsDiffs,1) ./ ...
                (nanstd(curRsDiffs,1) ./ sqrt(curSubjWEnoughVox));
            
            % Now take real tvals and null distribution of tvals and compare to get
            % pvals.
            pvalsDiffsEachRoiAndModelCombo(iRoi, iModelCombo) = getTwoTailedPermPval( ...
                tvalsDiffsEachRoiAndModelCombo(iRoi,iModelCombo), ...
                nullTvalsDiffsEachRoiAndModelCombo(iRoi,iModelCombo,:));
        end
    end
    
    % Now make/save pvals and stats
    % -----------------------------------------------------
    
    pvalsByModelAndRoi = array2table(pvalsEachRoiAndModel, ...
        'VariableNames',featLabels, 'RowNames', roiNamesToUse);
    disp(pvalsByModelAndRoi);
    
    modelComparisonPvals = array2table(pvalsDiffsEachRoiAndModelCombo, ...
        'VariableNames', genvarname(modelComboStrings), ...
        'RowNames', roiNamesToUse);
    disp(modelComparisonPvals);
end


%% Now also make summary plot across subjects/ROIs (Figure 3a).

if useSavedIntermedFiles && doMainAnalysis
    
    featBarColors = [ ...
        .7 .7 .2; ...
        .8 .5 .2; ...
        .7 .2 .6];
    featScatterColors = [ ...
        .5 .5 .1; ...
        .6 .4 .1; ...
        .5 .1 .4];
    voiFontSize = 12;
    axLabelFontSize = 12;
    
    if ~exist('rToPlotEachRoiModelAndSubj','var') || ~exist('nVoxBySubjAndRoi','var')
        error('Woops! You need to run the previous section first.')
    end
    plotThingsOnTopOfGroupedBarPlot(rToPlotEachRoiModelAndSubj, 3, 1, 2, ...
        pvalsEachRoiAndModel, featBarColors, featScatterColors);
    title({'Figure 3A','Model performance', 'up to 200 best localizer voxels','min 75'});
    ylabel('Mean of corrs w held-out data');
    % xlabel('ROI');
    xticklabels(roiNamesToUse); % , 'fontsize', voiFontSize);
    curAx = gca;
    curAx.XAxis.FontSize = voiFontSize;
    curAx.XAxisLocation = 'origin'; curAx.TickLength = [0 0];
    % curAx.XAxis.FontWeight = 'bold'; curAx.XLabel.FontWeight = 'normal';
    curAx.XLabel.FontSize = axLabelFontSize;
    curAx.YAxis.FontSize = axLabelFontSize; % (this does label and ticks this size)
    
    
    legend(featLabels, 'FontSize', axLabelFontSize);
    
end


%% Helpful functions that're ok to live in this script (i.e., don't need to be used across scripts)
%{
function curPrctile975 = plotMeanOnNullDistHist(meanRsByVox,meanRsByVox_perms,roiColor,titleText)
    figure; hold on;
    
    curNullVoxMeanArray = mean(meanRsByVox_perms,1); % Mean across rows, which are vox; leaves 10,000 null-dist means.
    curMean = mean(meanRsByVox);
    permsBelowMean = sum(curNullVoxMeanArray < curMean);
    % Commented lines are included in function now
    % lessExtremePerms = sum(abs(curNullVoxMeanArray) < abs(curMean));
    nPerms = size(meanRsByVox_perms,2);
    % curPval = 1 - ((lessExtremePerms-1)/ nPerms);
    curPval = getTwoTailedPermPval(curMean,curNullVoxMeanArray);
    curPrctile975 = prctile(curNullVoxMeanArray,97.5);
    % prctile975BySubjAndRoi(iSubj,iRoi) = curPrctile975;
    
    histogram(curNullVoxMeanArray,'FaceColor',[0.3 0.3 0.3], 'EdgeColor', [.4 .4 .4]);
    line([curMean curMean], ylim, 'Color', roiColor, ...
        'LineWidth', 3, 'LineStyle', ':');
    line([curPrctile975 curPrctile975], ylim, 'Color', [0.3 0.3 0.3], ...
        'LineStyle',':');
    curYLim = ylim;
    text(curMean - (curMean/50), curYLim(2)-(curYLim(2)/25), ...
        sprintf('%d/%d permutations\n below actual mean;\nTwo-tailed p < %0.4f', ...
        permsBelowMean, nPerms, curPval), 'FontWeight','bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment','top');
    xlabel('Prediction perf. (r), mean across vox');
    title(titleText);
    
    end
%}

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

function [sortedValsToCheck, voxInds] = getTopNVoxInds(curValsToRankBy, nVox)

[sortedValsToCheck, highestValInds] = ...
    sort(curValsToRankBy, 'descend');
voxInds = highestValInds(1:nVox);
end


function [curArrayOfAvgs, plotHandle] = plotVoxSelectionPlotWithinLoop(rsBySubjAndRoi, ...
    sortedCorrespondingTvals, highestCorrespondingTvalInds, ... % these two inputs should be controlled outside loop (either evc or scene tvals depending on ROI).
    iRoi, iSubj, curRoi, curSubjNum, plotCtr)

locLabelsByRoi = {'evc', 'evc', 'evc', 'scene', 'scene', 'scene'};


% Figure showing effects of vox selection
% -------------------------------------------------------

subplot(2,3,plotCtr); hold on;

curArrayOfAvgs = nan(numel(sortedCorrespondingTvals),1); % Size of number of voxels bc starting by including one, ending by including all

for iNHighest = 1:numel(sortedCorrespondingTvals)
    voxToInclude = highestCorrespondingTvalInds(1:iNHighest);
    curArrayOfAvgs(iNHighest) = ...
        mean(rsBySubjAndRoi{iSubj,iRoi}(voxToInclude));
end
plotHandle = plot(curArrayOfAvgs);
xlabel(['Top N ' locLabelsByRoi{iRoi} ' loc. voxels included in avg. perf.']);


ylabel('Prediction performance (r)');
curYlim = ylim;
line([75 75], ylim, 'Color', [0.8 0.8 0.8], 'LineStyle', ':');
line([200 200], ylim, 'Color', [0.3 0.3 0.3]); % max 200 vox used, even though more in arrays.
xlim([-0.5 200.5]); % 4.5]);
title(sprintf('%s, %s', curSubjNum, curRoi));

end
%}
