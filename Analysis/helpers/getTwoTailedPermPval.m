function curPval = getTwoTailedPermPval(curMean,nullDist)

% Get two-tailed p-value, without the assumption of a symmetrical null
% distribution (sees what proportion of the null values the actual mean is
% *more extreme* than)

lessExtremePerms = sum(abs(nullDist) < abs(curMean));
nPerms = numel(nullDist);
curPval = 1 - ((lessExtremePerms-1)/ nPerms);
