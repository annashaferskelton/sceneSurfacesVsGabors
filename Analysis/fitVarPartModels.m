%% This script fits models for variance partitioning (For figures 3B and 3C)
% Then will run visualization script after

%% Setup
clear;
close all;
% numCores = 8; % Max number of cores to use on the server.
% doGpu = 0; % Looks like it's much faster _not_ to gpu things, or else I'm doing it wrong!
% singleArrays = 0;
% excludeDepthArtifactIms = 0;
% veryFewLambdas = 1;
% doPermTests = 1;
nPerms = 100; % nPerms is the number of rand perms in addition to the actual analysis.
shorterTrWindow = 1;
% longerTrWindow = 0; % This stays this way
doAvgSignalAnalysis = 1;
% doFittingByTask = false;
% onlyLaterTrs = 0;
numCvSplitsOuter = 9; % 18; % 9; Can be 18 (2 runs grouped together) or 9 (4 runs grouped together); there are 36 runs total.
numCvSplitsInner = 10; % 8; % 10; % Was 10, but changing to 8 to keep same across-run structure as outer folds.

analysisDir = pwd;
filesepinds = find(analysisDir==filesep);
mainDir = analysisDir(1:filesepinds(end)-1);
behavDataDir = fullfile(mainDir,'DataBehavior');

subjIds = {'CI','CO','CK','CR','BX','CS','CG','CH'};
subjNums = {'S01','S02','S03','S04','S05','S06','S07','S08'};
subjsToDoNow = [1 2 3 4 5 6 7 8];

% thresholdedVois = 1; % This is the official option bc means the voxels passed a localizer threshold in independent localizer data.
% gssParcSceneAreas = 0;
%{
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

% arrayPrefix = '_double_';
%{
if longerTrWindow && shorterTrWindow
   error('Oh no! You can''t run script for multiple TR windows at once!');
end
if longerTrWindow % Controls which main signal file gets read as well as the ridge output file name
    avgTRs_targ = [6,10];
    trWindowText = '4-8sec';
elseif shorterTrWindow
    avgTRs_targ = [4,6];
%}
trWindowText = '2.4-4.8sec';
%{
else
    avgTRs_targ = [6,8]; % Will end up including 4-6.4 seconds after stim; includes times from future trials, as in Lescroart & Gallant (2019).
    trWindowText = '4-6ishSec';
end

if onlyLaterTrs
    trToStart = 6;
else
    trToStart = 1;
end
%}
if numCvSplitsOuter == 9 && numCvSplitsInner == 10
    cvFoldText = 'cv9';
elseif numCvSplitsOuter == 18 && numCvSplitsInner == 10
    cvFoldText = 'cv18';
elseif numCvSplitsOuter == 9 && numCvSplitsInner == 8
    cvFoldText = 'cv9-8';
end

% artifactPrefix = ''; % No artifact ims in this version bc ground-truth-ish depth
runsInSeparateCvFolds = 1; % <- this almost has to be the case.

%{
if veryFewLambdas
    lambdaText = '_fewerLambdas';
else
    lambdaText = '';
end
if doPermTests
    permText = '_permTests';
else
    permText = '';
end
%}

featLabels = {'spat3d','nonSpat3d','gabor'}; % also a pix model, but probably not involved in var part analysis.
featTypesA = [1, 2];
featTypesB = [3, 3];

for iFeatCombo = 1:numel(featTypesA)
    featTypeA = featTypesA(iFeatCombo);
    featTypeB = featTypesB(iFeatCombo);
    featPrefixA = ['_' featLabels{featTypeA} '_']; % e.g., '_nonSpat3d_'
    featDirA = fullfile(analysisDir,'ridgeInputs', ...
        [featLabels{featTypeA} 'Feats']);
    
    featPrefixB = ['_' featLabels{featTypeB} '_']; % e.g., '_gabor_'
    featDirB = fullfile(analysisDir,'ridgeInputs', ...
        [featLabels{featTypeB} 'Feats']);
    
    % outDirAB = fullfile(analysisDir, 'ridgeOutputs', ...
    % ['varPart' featPrefixA featPrefixB artifactPrefix lambdaText permText]);
    outDirAB = fullfile(analysisDir, 'ridgeOutputs', ...
        ['varPart' featPrefixA featLabels{featTypeB}]);
    if ~exist(outDirAB, 'dir')
        mkdir(outDirAB);
    end
    % Individual files will get saved depending on what portion of the loop
    % they're in! So saved below.
    
    %% Get permutations for doing permutation tests
    
    %% Now get features for each timepoint.
    tic;
    for iSubj = 1:numel(subjIds)
        curSubj = subjIds{iSubj};
        curSubjNumStr = subjNums{iSubj};
        % ----------------------------------------------------------------------
        % First, read in features for each trial
        % ----------------------------------------------------------------------
        
        % Get features for feat-types A and B
        
        curFeatFileA = fullfile(featDirA,['orderedFeats' featPrefixA curSubj '.mat']);
        load(curFeatFileA);
        
        if featTypeA == 2
            curFeatsA = curNonSpatFeats;
            clear curNonSpatFeats;
        elseif featTypeA == 1
            curFeatsA = curSpatFeats;
            clear spatFeats;
        end
        
        curFeatFileB = fullfile(featDirB,['orderedFeats' featPrefixB curSubj '.mat']);
        load(curFeatFileB);
        
        if featTypeB == 2
            curFeatsB = curNonSpatFeats;
            clear curNonSpatFeats;
        elseif featTypeB == 3
            curFeatsB = curFeats;
            clear curFeats;
        end
        
        % No exclusions here
        % curExclusionFile = fullfile('ridgeInputs',['orderedGoodTrialInds_' curSubj '.csv']);
        % tempGoodTrials = readtable(curExclusionFile,'ReadRowNames',false);
        % curGoodTrials = table2array(tempGoodTrials(:,2));
        % imNamesToCompare2 = tempGoodTrials.Var1;
        
        % if ~isequal(imNamesToCompare,imNamesToCompare2); error('oh no! Im orders don''t match!'); end
        %{
    if runsInSeparateCvFolds
        % Read in the run info here!!!
        trialDatFile = fullfile('ridgeInputs',['orderedScenesOnlyTrialDat_' curSubj '.csv']);
        tempTrialDat = readtable(trialDatFile,'ReadRowNames',false);
        
        imNamesToCompare3 = tempTrialDat.ImgName;
        % if ~isequal(imNamesToCompare,imNamesToCompare3); error('oh no! Im orders don''t match!'); end
        % spot-checked orders and looks good, but I think .jpg and .jpeg are
        % messing things up.
        
        if excludeDepthArtifactIms ==1
            trialDat = tempTrialDat(curGoodTrials==1,:);
        else
            trialDat = tempTrialDat;
        end
        sessInfo = trialDat.Sess;
    else
        sessInfo = [];
    end
        %}
        
        % Exclude from z-scoring features in this case.
        %{
    if excludeDepthArtifactIms == 1
        curFeatsA = curFeatsA(curGoodTrials==1,:);
        curFeatsB = curFeatsB(curGoodTrials==1,:);
    end
        %}
        
        % Read in the run info here!!!
        trialDatFile = fullfile(behavDataDir,curSubjNumStr,'intermedArrays', ...
            [curSubjNumStr '_byRunAndTrial.mat']);
        load(trialDatFile,'byRun','byRunAndTrial');
        runNums = reshape(byRunAndTrial.globalRunNum',[],1);
        taskNums = reshape(byRunAndTrial.task',[],1);
        
        
        % Z-score features
        % --------------------------------------------------------
        % ("all feature channels were normalized to have a mean of
        % 0 and a standard deviation of 1 (z-scored)")
        curFeatsA_zd = zscore(curFeatsA, [], 1); % becomes feats by trials; gabor feats should already be z-scored
        figure; imagesc(curFeatsA_zd); title(sprintf('%s feats',featPrefixA));
        figure; histogram(curFeatsA_zd); title(sprintf('distribution of %s feats',featPrefixA));
        designMatA = curFeatsA_zd;
        
        curFeatsB_zd = zscore(curFeatsB, [], 1); % becomes feats by trials; gabor feats should already be z-scored
        figure; imagesc(curFeatsB_zd); title(sprintf('%s feats',featPrefixB));
        figure; histogram(curFeatsB_zd); title(sprintf('distribution of %s feats',featPrefixB));
        designMatB = curFeatsB_zd;
        
        fprintf('Done processing features for %s!\n', curSubj);
        
        %% Set up for model selection within model selection data.
        
        % Make 9 blocks of 4 runs each OR 18 blocks of 2 runs each to
        % keep together when running ridge code -- original option set
        % at top of script!
        nRuns = numel(unique(runNums));
        nRunsGrouped = nRuns / numCvSplitsOuter;
        groupedRuns = ceil(runNums/nRunsGrouped); % groupedRuns array will be the same also for task breakdowns, just subsetted every other run (each grouping half as big)
        % ^ Groups to use for cross-validation
        
        % lambdas = [0 logspace(-100,10)]; % for now (so far, subj 2), lambdas log-spaced. Check if good decision!
        % if veryFewLambdas
        % lambdas = 0;
        lambdas = [0 logspace(-2,5,12)];
        % else
        % lambdas = logspace(-15,7,200);% 0:0.1:1000;
        % end
        
        
        
        %% Now read in the data I need for this analysis.
        
        % curDataFile = fullfile(mainDir,'Samples',['MainTaskSignalByTrial_' trWindowText sceneText '_', curSubjNumStr '.mat']);
        curDataFile = fullfile(mainDir,'DataFmri',['MainTaskSignalByTrial_' trWindowText '_', curSubjNumStr '.mat']);
        
        load(curDataFile,'ROI_names','mainSig');
        endRoi = numel(ROI_names)
        
        
        for iRoi = roisToDoNow
            
            curRoi = roiNames{iRoi};
            
            if doAvgSignalAnalysis % && ~doFittingByTask
                % Read in ridge inputs!
                % -------------------------------------
                
                curSubjDat = mainSig(iRoi).dat_avg';
                
                figure; imagesc(curSubjDat);
                disp('Done reading in and displaying data!');
                
                % perm array should be perms x trials x outer splits
                % set seed so same order for subjects/ROIs (and doing it out here
                % bc using same for models and concatenated feats)
                rng('default');
                randSeedArray = nan(nPerms,numCvSplitsOuter);
                numVoxToModel = size(curSubjDat,1);
                counter = 0;
                for iSplit=1:numCvSplitsOuter
                    for iPerm = 1:nPerms
                        counter = counter+1;
                        randSeedArray(iPerm,iSplit) = counter;
                    end
                end
                
                % -----------
                % Start ridge code
                % -------
                % perm array should be perms x trials x outer splits
                % Separate file saved for every subj/roi/model (perm array same across subjs/rois).
                % permArray = randperm();
                % tic;
                %{
            [curPredRsA,curPredRsA_perms, curBsA,lambdasBySplitA] = ...
                doRidge_cvLambdas_wRandPermsInput_noMean(curSubjDat,designMatA,numCvSplitsOuter,numCvSplitsInner,lambdas, ...
                doGpu,numCores,singleArrays,runsInSeparateCvFolds, randSeedArray, groupedRuns); %,numCores);
            [curPredRsB,curPredRsB_perms, curBsB,lambdasBySplitB] = ...
                doRidge_cvLambdas_wRandPermsInput_noMean(curSubjDat,designMatB,numCvSplitsOuter,numCvSplitsInner,lambdas, ...
                doGpu,numCores,singleArrays,runsInSeparateCvFolds, randSeedArray, groupedRuns);
            
            designMatConcat = [designMatA designMatB];
            [curPredRsConcat,curPredRsConcat_perms, curBsConcat,lambdasBySplitConcat] = ...
                doRidge_cvLambdas_wRandPermsInput_noMean(curSubjDat,designMatConcat,numCvSplitsOuter,numCvSplitsInner,lambdas, ...
                doGpu,numCores,singleArrays,runsInSeparateCvFolds, randSeedArray, groupedRuns);
                %}
                
                [curPredRsA,curPredRsA_perms, curBsA,lambdasBySplitA] = ...
                    doRidge_cvLambdas_wRandPermsInput_noMean(curSubjDat,designMatA,numCvSplitsOuter,numCvSplitsInner,lambdas, ...
                    runsInSeparateCvFolds, randSeedArray, groupedRuns); %,numCores);
                [curPredRsB,curPredRsB_perms, curBsB,lambdasBySplitB] = ...
                    doRidge_cvLambdas_wRandPermsInput_noMean(curSubjDat,designMatB,numCvSplitsOuter,numCvSplitsInner,lambdas, ...
                    runsInSeparateCvFolds, randSeedArray, groupedRuns);
                
                designMatConcat = [designMatA designMatB];
                [curPredRsConcat,curPredRsConcat_perms, curBsConcat,lambdasBySplitConcat] = ...
                    doRidge_cvLambdas_wRandPermsInput_noMean(curSubjDat,designMatConcat,numCvSplitsOuter,numCvSplitsInner,lambdas, ...
                    runsInSeparateCvFolds, randSeedArray, groupedRuns);
                
                % toc;
                %{
            fprintf('Done fitting variance partitioning %s,%s models, sub %s, roi %s!\n', featPrefixA, featPrefixB, curSubj, curRoi);
            save(fullfile(outDirAB,['fittedModelsConcat' featPrefixA featPrefixB curSubj '_' curRoi sceneText arrayPrefix '_vectorized_collapsedHemispheres.mat']));
                %}
                
                fprintf('Done fitting variance partitioning %s,%s models, sub %s, roi %s!\n', featPrefixA, featPrefixB, curSubj, curRoi);
                save(fullfile(outDirAB,['fittedModelsConcat' featPrefixA featPrefixB curSubj '_' curRoi '_collapsedHemispheres.mat']));
                
                % First, histogram of the mean r's by voxel.
                % figure; histogram(meanRsByVoxA); xlabel('r''s'); title(sprintf('Mean correlations, %s, %s %s models, %s voxels, %s', artifactPrefix(2:end), featPrefixA, featPrefixB, curRoiCollapsed, curSubj));
                
            end
        end; clear curRoi; % iRoi
    end; clear curSubj curSubjNumStr; % iSubj
end
toc;
%% diagnostic plots are now in separate script


%% Separate ridge func
function [curPredRs, curPredRs_perms, curBs,lambdasBySplitAndVox] = ...
    doRidge_cvLambdas_wRandPermsInput_noMean(curSubjDat,designMat,numCvSplitsOuter,numCvSplitsInner, ...
    lambdas, doCvByRuns, randSeedArray, varargin) %, numCores) % Last arg is run info!

doGpu = 0;
isSingle = 0;

% perm array should be perms x trials x outer splits
nPerms = size(randSeedArray,1); % rows in perm array are a list of shuffled labels. (as many rows as perms, as many cols as trials)

if doCvByRuns == 1
    eachTrialsRun = varargin{1}; % This is what's used to group into cv folds, regardless of which is actually a "run"
elseif doCvByRuns == 0
else
    error('oh no! doCvByRuns has to either be 1 or 0.');
end

numVox = size(curSubjDat,1);
numTrsTrn = size(curSubjDat,2);
% numTrsHeldOut = numTrsTrn / numCvSplits; % For now, TRs = Trials! and ignoring fact that it says 'trn'
numPreds = size(designMat,2); % + 1; %  <- this is defined inside the function as well as outside.

lambdasBySplitAndVox = nan(numCvSplitsOuter,numVox);
curBs = nan(numPreds,numVox,numCvSplitsOuter);
curPredRs = nan(numVox,numCvSplitsOuter);
curPredRs_perms = nan(numVox,numCvSplitsOuter,nPerms);

if isSingle
    % Should only need to
    lambdasBySplitAndVox = single(lambdasBySplitAndVox);
    curBs = single(curBs);
    curPredRs = single(curPredRs);
    curSubjDat = single(curSubjDat);
    designMat = single(designMat);
    lambdas = single(lambdas);
end
if doGpu
    lambdasBySplitAndVox = gpuArray(lambdasBySplitAndVox);
    curBs = gpuArray(curBs);
    curPredRs = gpuArray(curPredRs);
    curSubjDat = gpuArray(curSubjDat);
    designMat = gpuArray(designMat);
end

for curSplit = 1:numCvSplitsOuter % For comparing performance of models (that we already think have the optimal k's chosen, as measured by predicting within each of these folds.)
    
    fprintf('Starting cross-validation split #%d \n', curSplit);
    
    if doCvByRuns == 1 % Keeping different sessions (and therefore runs) separate across cv folds
        isHeldOut = cvPartitionsByRunCurSplit(numCvSplitsOuter,eachTrialsRun,curSplit);
    else % Evenly divided cv
        isHeldOut = evenCvPartitionsCurSplit(numTrsTrn,numCvSplitsOuter,curSplit);
    end
    
    % Separate into held-out and not for the data and also features
    % from each model. (for main CV loop)
    curY = curSubjDat(:,~isHeldOut)'; % Transposing so it's a column vector. (Now should be TRs by vox.)
    curHeldOutY = curSubjDat(:,isHeldOut==1)'; % Transposing so it's a column vector. (Now should be all vox.)
    
    % Get predictors for used vs. held-out data on this split.
    curX = designMat(~isHeldOut,:); % [ones(sum(~isHeldOut),1) designMat(~isHeldOut,:)]; % transpose bc observations are rows, features are columns. Then put column for constant.
    curHeldOutX = designMat(isHeldOut==1,:); % [ones(sum(isHeldOut),1) designMat(isHeldOut==1,:)];
    
    
    % TO STORE THINGS FROM INNER CV LOOP (use to pick best lambdas
    % during rest of this larger cv fold)
    % -----------------------------------------------------------------
    
    curInnerRSqs = ...
        nan(numel(lambdas),numVox,numCvSplitsInner); % To store r-squared between predicted and actual held-out data of each *inner* cv loop.
    if isSingle
        curInnerRSqs = single(curInnerRSqs);
    end
    if doGpu
        curInnerRSqs = gpuArray(curInnerRSqs);
    end
    
    % Inner split of the data, used to choose lambdas that generalize, then comparing best-case model fits for each model type.
    for curSplitInner = 1:numCvSplitsInner
        
        % For indexing into 9/10 of this split's data.
        numTrsInner = sum(~isHeldOut);
        
        if mod(numTrsInner,numCvSplitsInner) > 0
            isHeldOutInner = almostEvenCvPartitionsCurSplit(numTrsInner,numCvSplitsInner,curSplitInner);
        elseif mod(numTrsInner,numCvSplitsInner) == 0
            isHeldOutInner = evenCvPartitionsCurSplit(numTrsInner,numCvSplitsInner,curSplitInner);
        end
        
        % Separate into held-out and not for the data and also features
        % from each model. (for inner CV loop)
        curYInner = curY(~isHeldOutInner,:); % In both cases, index into the not-held-out data from the outer loop.
        curHeldOutYInner = curY(isHeldOutInner==1,:);
        
        % Get predictors for used vs. held-out data on this split.
        % Always taking from not-held-out data from the outer loop.
        curXInner = curX(~isHeldOutInner,:);
        curHeldOutXInner = curX(isHeldOutInner==1,:);
        
        
        % parfor (iLambda = 1:numel(lambdas),numCores) % Max of 8 cores on server (to be nice; can change if nobody using)
        % parfor (iLambda = 1:numel(lambdas),numCores)
        for iLambda = 1:numel(lambdas)
            curLambda = lambdas(iLambda);
            
            % tic;
            curBsInner = ...
                (curXInner'*curXInner + curLambda*eye(numPreds)) \ curXInner'*curYInner;
            
            curYHatsInner = curHeldOutXInner * curBsInner;
            curCorrInner = diag(corr(curYHatsInner,curHeldOutYInner)); % r-squared between predicted and held-out data
            curInnerRSqs(iLambda,:,curSplitInner) = (curCorrInner .^ 2) .* sign(curCorrInner); % Get square for each vox, save.
            
        end % iLambda
        
    end % curSplitInner
    
    if rand < 0.1
        plotInnerRSqs(curInnerRSqs,lambdas,15,'10 perc trn'); % 15 is how many vox
    end
    
    % Find best lambda (all voxels together)
    % -----------------------------------------------------
    curAccsByLambdaVox = mean(curInnerRSqs,3); % mean of r-squareds over all the inner splits. lambda by vox
    % Store most successful weights for this voxel.
    [~,bestLambdaInds] = max(curAccsByLambdaVox,[],1); % If 2 are equally good, it'll take the smaller lambda, which we want. (should be over dim 1, for each vox.)
    bestLambdas = lambdas(bestLambdaInds);
    lambdasBySplitAndVox(curSplit,:) = bestLambdas; % Best lambdas for this outer split
    
    % Now find final weights for best lambda and this split (since each vox has a diff lambda, I think I need a loop!)
    
    for iVox = 1:numVox % (doing separately for each vox here)
        curBs(:,iVox,curSplit) = ...
            (curX'*curX + bestLambdas(iVox)*eye(numPreds)) \ curX'*curY(:,iVox);
    end
    curYHats = curHeldOutX * curBs(:,:,curSplit);
    curPredRs(:,curSplit) = diag(corr(curYHats,curHeldOutY)); % Pearson's r between predicted and held-out data
    
    if nPerms > 0
        for iPerm = 1:nPerms
            % Shuffle the different predicted image features... (verify
            % later this makes the most sense???)
            rng(randSeedArray(iPerm,curSplit)); % Set to rand seed that's the same for this perm, this outer split, for all subjs/rois
            permYHats = curYHats(randperm(size(curYHats,1)),:); % Each column is a voxel; randomly shuffle the rows.
            curPredRs_perms(:,curSplit,iPerm) = diag(corr(permYHats,curHeldOutY));
        end
    end
    
end % curSplit
% meanRsByVox = mean(curPredRs,2);
% meanRsByVox_perms = squeeze(mean(curPredRs_perms,2));
end

%% Define function for cv partitions

function isHeldOut = evenCvPartitionsCurSplit(nTrs, nSplits, curSplit)
nTrsHeldOut = nTrs / nSplits; % Split 90% of data (not-held-out) into 10 pieces.

isHeldOut = zeros(1,nTrs);
heldOutTrsInnerStartInd = ((curSplit-1) * nTrsHeldOut) + 1;
heldOutTrsInnerEndInd = curSplit * nTrsHeldOut;
isHeldOut(heldOutTrsInnerStartInd:heldOutTrsInnerEndInd) = 1;
end

function isHeldOut = almostEvenCvPartitionsCurSplit(nTrs, nSplits, curSplit)

nTrsHeldOutMost = floor(nTrs / nSplits); % Split 90% of data (not-held-out) into 10 pieces.
% nTrsHeldOutLast = nTrs - (nTrsHeldOutMost*(nSplits-1));
isHeldOut = zeros(1,nTrs);
heldOutTrsStartInd = ((curSplit-1) * nTrsHeldOutMost) + 1; % Start point doesn't change depending on cv fold (since last fold is odd one out)

if curSplit < nSplits
    heldOutTrsEndInd = curSplit * nTrsHeldOutMost;
elseif curSplit == nSplits
    heldOutTrsEndInd = nTrs;
end

isHeldOut(heldOutTrsStartInd:heldOutTrsEndInd) = 1;
disp(['CV fold ' num2str(curSplit) '/' num2str(nSplits) ', ' num2str(sum(isHeldOut)) '/' num2str(length(isHeldOut)) ' held out TRs'])

end


function isHeldOut = cvPartitionsByRunCurSplit(nSplits, eachTrialsRun, curSplit)
nParts = length(unique(eachTrialsRun));
nPartsHeldOut = nParts / nSplits; % Split 90% of data (not-held-out) into 10 pieces.

heldOutPartStartInd = ((curSplit-1) * nPartsHeldOut) + 1;
heldOutPartEndInd = curSplit * nPartsHeldOut;
curHeldOutParts = heldOutPartStartInd:heldOutPartEndInd;

isHeldOut = ismember(eachTrialsRun,curHeldOutParts);
disp(['CV fold ' num2str(curSplit) '/' num2str(nSplits) ', ' num2str(sum(isHeldOut)) '/' num2str(length(isHeldOut)) ' held out TRs'])

end

function plotInnerRSqs(rSqLambdaVoxFold,lambdas,nVox,cvString)
s = RandStream('mlfg6331_64');
randVox = randsample(s,1:size(rSqLambdaVoxFold,2), nVox); % nVox is how many to select
sAcross = ceil(sqrt(nVox));
sDown = ceil(nVox/sAcross);
figure;
for iVox = 1:nVox
    curVox = randVox(iVox);
    subplot(sDown,sAcross,iVox);
    % plot(lambdas,squeeze(rSqLambdaVoxFold(:,curVox,:)));
    semilogx(lambdas,squeeze(rSqLambdaVoxFold(:,curVox,:)),'Color',[0.3 0.3 0.3]);
    
    hold on; box off;
    bestLambdas = [];
    for iFold = 1:size(rSqLambdaVoxFold,3)
        [~,maxInd] = max(squeeze(rSqLambdaVoxFold(:,curVox,iFold)));
        bestLambdas(iFold) = lambdas(maxInd);
        % curMaxLine = line([lambdas(maxInd),lambdas(maxInd)],ylim,'LineWidth',2,'Color',[0.3 0.3 0.3], 'LineStyle',':');
        % curMaxLine.Color(4) = 0.3;
    end
    
    
    p = semilogx(lambdas,squeeze(mean(rSqLambdaVoxFold(:,curVox,:),3)), ...
        'LineWidth',2,'Color','red');
    p.Color(4) = 0.5;
    
    [~,maxInd] = max(squeeze(mean(rSqLambdaVoxFold(:,curVox,:),3)));
    pAvg = line([lambdas(maxInd),lambdas(maxInd)],ylim, 'LineWidth',2,'Color','red','LineStyle',':');
    pAvg.Color(4) = 0.5;
    
    % Also get lambda with highest accurace *across* splits of data.
    maxEachFold = max(squeeze(rSqLambdaVoxFold(:,curVox,:))); % Max over first dimension (lambda)
    [curMax,bestFold] = max(maxEachFold);
    [curMaxAgain,bestLBestFold] = max(rSqLambdaVoxFold(:,curVox,bestFold));
    o = line([lambdas(bestLBestFold),lambdas(bestLBestFold)],ylim,'LineWidth',2,'Color','magenta','LineStyle',':');
    o.Color(4) = 0.7;
    
    title(sprintf('Vox %d, %s',curVox,cvString));
    %}
    
end
sgtitle('Best lambdas across folds');
end
