function CommonResponseTest = frstat_common_response_test( ...
    est_all, y_grid, evecs, evals, evecsx, evalsx, XX2, sampleSize, ...
    bmgrid, nRep, seed, betaDraws)
%FRSTAT_COMMON_RESPONSE_TEST Test equality of the two CLR response functions.
%
% The null is beta_1(s) = beta_2(s) in the common physical forcing unit
% (W/m^2). The statistic is the effective sample size times the L2 norm of
% the estimated coefficient-function difference. Critical values use the plug-in limit simulation
% as the forcing-specific response inference, applied directly to e_1-e_2.

if size(est_all, 2) ~= 2 || size(XX2, 2) ~= 2
    error('The common-response test requires exactly two forcing coordinates.');
end
if nRep < 1 || bmgrid < 2
    error('nRep must be positive and bmgrid must be at least two.');
end

y_grid = y_grid(:);
contrast = [1; -1];
estimatedDifference = est_all(:, 1) - est_all(:, 2);
observedStatistic = sampleSize * sqrt(trapz(y_grid, estimatedDifference.^2));

if nargin < 12 || isempty(betaDraws)
    betaDraws = frstat_simulate_beta_draws(2, numel(evals), bmgrid, XX2, ...
        evecs, evals, evecsx, evalsx, sampleSize, nRep, seed);
else
    if size(betaDraws, 1) ~= numel(y_grid) || size(betaDraws, 2) ~= 2
        error('betaDraws must have dimensions length(y_grid)-by-2-by-nRep.');
    end
    nRep = size(betaDraws, 3);
end

contrastDraws = sampleSize * reshape(betaDraws(:,1,:) - betaDraws(:,2,:), ...
    numel(y_grid), nRep);
simulatedStatistics = sqrt(trapz(y_grid, contrastDraws.^2, 1))';

valid = isfinite(simulatedStatistics);
if nnz(valid) ~= nRep
    error('frstat_common_response_test:InvalidMonteCarloDraws', ...
        'Only %d of %d requested coefficient draws were finite.', ...
        nnz(valid), nRep);
end
simulatedStatistics = simulatedStatistics(valid);

CommonResponseTest = struct();
CommonResponseTest.null = "beta_F1(s) = beta_F2(s) in W/m^2 units";
CommonResponseTest.contrast = contrast;
CommonResponseTest.observedStatistic = observedStatistic;
CommonResponseTest.criticalValue90 = quantile(simulatedStatistics, 0.90);
CommonResponseTest.criticalValue95 = quantile(simulatedStatistics, 0.95);
CommonResponseTest.criticalValue99 = quantile(simulatedStatistics, 0.99);
CommonResponseTest.pValue = ...
    (1 + sum(simulatedStatistics >= observedStatistic)) / ...
    (numel(simulatedStatistics) + 1);
CommonResponseTest.nValid = numel(simulatedStatistics);
CommonResponseTest.seed = seed;
CommonResponseTest.bmgrid = bmgrid;
CommonResponseTest.estimatedDifference = estimatedDifference;
end
