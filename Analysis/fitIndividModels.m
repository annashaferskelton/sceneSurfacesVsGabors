%% This version of code separates ridge fitting into a separate function
% Also makes sure I've streamlined stuff (e.g., shape of design matrices)
% to keep from switching things then switching them back

% 2024.05.31 cleaning up code for public repository

%% Setup
clear;
close all;
nPerms = 100; % 10000 was the number in the manuscript % nPerms is the number of rand perms in addition to the non-shuffled analysis.

doMainAnalysis = 1; % Analysis averaging across peak TRs instead of going TR by TR
numCvSplitsOuter = 9; % There are 36 runs total (12 in each of 3 sessions), so 9 cv splits will create groups of 4 consecutive runs
numCvSplitsInner = 10; % Used to choose ridge parameter within training data

% Analysis starts in this directory; find path of overall project directory
% without changing directories (to be able to grab behavioral data, e.g.,
% from its separate folder)
analysisDir = pwd;
filesepinds = find(analysisDir==filesep);
mainDir = analysisDir(1:filesepinds(end)-1);
behavDataDir = fullfile(mainDir,'DataBehavior');
fmriDataDir = fullfile(mainDir,'DataFmri');

% Participants have anonymized initials (subjIds) in addition to
% participant numbers (subjNums)
subjIds = {'CI','CO','CK','CR','BX','CS','CG','CH'};
subjNums = {'S01','S02','S03','S04','S05','S06','S07','S08'};
subjsToDoNow = [1:8]; % Indexes into above participant lists to only do the specified participant(s)


roiNames = {'V1','V2','V3', ...
    'OPA', 'PPA', 'RSC'};
roisToDoNow = 1:6; % Indexes into list above

% trWindowText tells you which TRs get included in this signal.
trWindowText = '2.4-4.8sec';

if numCvSplitsOuter == 9 && numCvSplitsInner == 10
    cvFoldText = 'cv9';
end

runsInSeparateCvFolds = 1; % <- this almost has to be the case.
featLabels = {'spat3d','nonSpat3d','gabor'};
featTypesToDoNow = 1:3; 

%% Now get features for each trial.
tic;
for iSubj = subjsToDoNow % 1:numel(subjIds)
    curSubj = subjIds{iSubj};
    curSubjNumStr = subjNums{iSubj};
    
    % ----------------------------------------------------------------------
    % First, read in features for each trial
    % ----------------------------------------------------------------------
    
    for featType = featTypesToDoNow 
        
        if featType == 1 % spatial 3d feats
            featPrefix = '_spat3d_';
            featDir = fullfile(analysisDir,'ridgeInputs','spat3dFeats');
            outDir = fullfile(analysisDir,'ridgeOutputs', ...
                ['spat3dFeats_' trWindowText '_' cvFoldText]);
            curFeatFile = fullfile(featDir,['orderedFeats' featPrefix curSubj '.mat']);
            load(curFeatFile,'curSpatFeats');
            curFeats = curSpatFeats;
        elseif featType == 2 % Non-spatial 3d feats
            featPrefix = '_nonSpat3d_';
            featDir = fullfile(analysisDir,'ridgeInputs','nonSpat3dFeats');
            outDir = fullfile(analysisDir,'ridgeOutputs', ...
                ['nonSpat3dFeats_' trWindowText '_' cvFoldText]);
            curFeatFile = fullfile(featDir,['orderedFeats' featPrefix curSubj '.mat']);
            load(curFeatFile,'curNonSpatFeats');
            curFeats = curNonSpatFeats;
        elseif featType ==3 % gabors
            featPrefix = '_gabor_';
            featDir = fullfile(analysisDir,'ridgeInputs','gaborFeats');
            outDir = fullfile(analysisDir,'ridgeOutputs', ...
                ['gaborFeats_' trWindowText '_' cvFoldText]); % fakeDataText]);
            curFeatFile = fullfile(featDir,['orderedFeats' featPrefix curSubj '.mat']);
            
            load(curFeatFile);
        else
            error('This feat type not done yet!');
        end
        
        if ~exist(outDir, 'dir')
            mkdir(outDir)
        end
        
        % Read in the run info here!!!
        trialDatFile = fullfile(behavDataDir,curSubjNumStr,'intermedArrays', ...
            [curSubjNumStr '_byRunAndTrial.mat']);
        load(trialDatFile,'byRun','byRunAndTrial');
        runNums = reshape(byRunAndTrial.globalRunNum',[],1);
        taskNums = reshape(byRunAndTrial.task',[],1);
        if numel(byRun.accs) < 36
            error('Oh no! Not all data''s here! Run behav analysis with all 3 sessions.')
        end
        
        % Z-score features
        % --------------------------------------------------------
        
        curFeats_zd = zscore(curFeats, [], 1); % becomes feats by trials; gabor feats should already be z-scored
        % figure; imagesc(curFeats_zd);
        designMat = curFeats_zd;
        
        fprintf('Done processing features for %s!\n', curSubj);
        
        %% Set up for model selection within model selection data.
        
        % Make 9 blocks of 4 runs each OR 18 blocks of 2 runs each to
        % keep together when running ridge code -- original option set
        % at top of script!
        nRuns = numel(unique(runNums));
        nRunsGrouped = nRuns / numCvSplitsOuter;
        groupedRuns = ceil(runNums/nRunsGrouped); % groupedRuns array will be the same also for task breakdowns, just subsetted every other run (each grouping half as big)
        
        lambdas = [0 logspace(-2,5,12)];
        
        %% Now read in the data I need for this analysis.
        
        curDataFile = fullfile(fmriDataDir,['MainTaskSignalByTrial_' trWindowText '_', curSubjNumStr '.mat']);
        load(curDataFile,'ROI_names','mainSig');
        endRoi = numel(ROI_names)
        
        for iRoi = roisToDoNow % 1:endRoi % numel(roiNames)
            
            curRoi = roiNames{iRoi};
            
            if doMainAnalysis 
                
                % Read in ridge inputs
                % -------------------------------------
                
                curSubjDat = mainSig(iRoi).dat_avg';
                
                disp('Done reading in and displaying data!');
                numVoxToModel = size(curSubjDat,1);
                
                % -----------
                % Start ridge code
                % -------
                
                [meanRsByVox,meanRsByVox_perms, curBs,lambdasBySplit] = ...
                    doRidge_cvLambdas_wRandPerms_noInt(curSubjDat,designMat,numCvSplitsOuter,numCvSplitsInner,lambdas, ...
                    runsInSeparateCvFolds, nPerms, groupedRuns); %,numCores);
                
                if isnan(mean(meanRsByVox))
                    keyboard;
                end
                
                fprintf('Done fitting %s models, sub %s, roi %s!\n', featPrefix, curSubj, curRoi);
                save(fullfile(outDir, ...
                    ['fittedModels' featPrefix curSubj '_' curRoi ... arrayPrefix ...
                    '_collapsedHemispheres_' trWindowText ...
                    '_' cvFoldText '.mat'])); % fakeDataText sceneText 
                
                % First, histogram of the mean r's by voxel.
                % figure; histogram(meanRsByVox); xlabel('r''s'); title(sprintf('Mean correlations, %s, %s models, %s voxels, %s', artifactPrefix(2:end), featPrefix, curRoi, curSubj));
            end
            
        end % iRoi
    end % featType
end % iSubj
toc;
%% diagnostic plots now in separate script


%% Separate ridge func
function [meanRsByVox, meanRsByVox_perms, curBs,lambdasBySplitAndVox] = ...
    doRidge_cvLambdas_wRandPerms_noInt(curSubjDat,designMat,numCvSplitsOuter,numCvSplitsInner, ...
    lambdas, doCvByRuns, nPerms, varargin) %, numCores) % Last arg is run info!

isSingle = 0;
doGpu = 0;
if doCvByRuns == 1
    eachTrialsRun = varargin{1};
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
        
        for iLambda = 1:numel(lambdas)
            curLambda = lambdas(iLambda);
            
            curBsInner = ...
                (curXInner'*curXInner + curLambda*eye(numPreds)) \ curXInner'*curYInner;
            
            curYHatsInner = curHeldOutXInner * curBsInner;
            curCorrInner = diag(corr(curYHatsInner,curHeldOutYInner)); % r-squared between predicted and held-out data
            curInnerRSqs(iLambda,:,curSplitInner) = (curCorrInner .^ 2) .* sign(curCorrInner); % Get square for each vox, save.
            fprintf('  outer cv fold %d/%d, inner split %d/%d, inner lambda %d/%d\n', ...
                curSplit,numCvSplitsOuter,curSplitInner,numCvSplitsInner, ...
                iLambda,numel(lambdas));
        end % iLambda
        
    end % curSplitInner
    
    
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
    
    % Make null distribution
    if nPerms > 0
        for iPerm = 1:nPerms
            % Shuffle the different predicted image features.
            permYHats = curYHats(randperm(size(curYHats,1)),:); % Each column is a voxel; randomly shuffle the rows.
            curPredRs_perms(:,curSplit,iPerm) = diag(corr(permYHats,curHeldOutY));
            fprintf('     Cur outer split %d/%d, cur perm %d/%d\n', ...
                curSplit,numCvSplitsOuter,iPerm,nPerms);
        end
    end
    
end % curSplit
meanRsByVox = mean(curPredRs,2);
meanRsByVox_perms = squeeze(mean(curPredRs_perms,2));
end

%% Define functions for cv partitions, other smaller tasks


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
    
end
sgtitle('Best lambdas across folds');
end
