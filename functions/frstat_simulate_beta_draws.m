function draws = frstat_simulate_beta_draws(nForce, KK, bmgrid, XX2, ...
    evecs, evals, evecsx, evalsx, T, nRep, seed)
%FRSTAT_SIMULATE_BETA_DRAWS Joint draws of coefficient-function estimation error.

if nargin >= 11 && ~isempty(seed)
    priorRng = rng;
    cleanupRng = onCleanup(@() rng(priorRng));
    rng(seed, 'twister');
end

if nForce ~= size(XX2, 2) || nForce ~= size(evecsx, 1)
    error('The requested forcing dimension does not match the predictor system.');
end
if KK ~= numel(evals) || KK ~= size(evecs, 2)
    error('KK, evals, and evecs have inconsistent response dimensions.');
end
if bmgrid < 2 || nRep < 1 || T < 2
    error('bmgrid, nRep, and T must define a nondegenerate simulation.');
end

sdUN = sqrt(max(evals(:), 0));
sdX  = sqrt(max(evalsx(:), 0));
tolUN = 1e-10 * max(1,max(abs(evals(:))));
tolX = 1e-10 * max(1,max(abs(evalsx(:))));
if min(evals(:)) < -tolUN || min(evalsx(:)) < -tolX
    error('Limit covariance eigenvalues contain a materially negative value.');
end
nGrid = size(evecs, 1);
nPredictors = size(XX2, 2);
draws = zeros(nGrid, nForce, nRep);

for rep = 1:nRep
    responsePaths = zeros(KK, bmgrid + 1);
    for k = 1:KK
        responsePaths(k, :) = brownian_motion(0, 1, bmgrid, 0);
    end

    motionPaths = zeros(nPredictors, bmgrid + 1);
    for k = 1:nPredictors
        motionPaths(k, :) = brownian_motion(0, 1, bmgrid, 0);
    end

    responseLimit = evecs * (sdUN .* responsePaths);
    predictorLimit = evecsx * (sdX .* motionPaths);

    dResponseLimit = diff(responseLimit, 1, 2);
    predictorLimit = predictorLimit(:, 1:end-1);
    predictorLimit = predictorLimit - mean(predictorLimit, 2);
    denominator = (predictorLimit * predictorLimit') / bmgrid;
    if rcond(denominator) < 1e-12
        error('A simulated predictor denominator is numerically singular.');
    end
    operatorLimit = (denominator \ (predictorLimit * dResponseLimit'))';

    draws(:, :, rep) = operatorLimit(:, 1:nForce) / T;
end
end
