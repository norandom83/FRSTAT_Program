function DensityCI = frstat_density_ci(y_grid, f_ref, est_all, ...
    errorDraws, shockDeltaX, alpha)

y_grid = y_grid(:);
f_ref = normalize_density_trapz(max(f_ref(:), eps), y_grid);
nGrid = numel(y_grid);
nForce = size(est_all, 2);

if size(est_all, 1) ~= nGrid || size(errorDraws, 1) ~= nGrid || ...
        size(errorDraws, 2) ~= nForce
    error('Density-interval inputs have inconsistent dimensions.');
end
shockDeltaX = shockDeltaX(:);
if numel(shockDeltaX) ~= nForce || any(~isfinite(shockDeltaX))
    error('One finite shock magnitude is required for each forcing.');
end
if ~isscalar(alpha) || alpha <= 0 || alpha >= 0.5
    error('alpha must be strictly between zero and one half.');
end

nRep = size(errorDraws, 3);
logReference = log(max(f_ref, 1e-12));
DensityCI = repmat(struct('estimate', [], 'lower', [], 'upper', [], ...
    'alpha', alpha, 'nRep', nRep), nForce, 1);

for ip = 1:nForce
    logEstimate = logReference + shockDeltaX(ip)*est_all(:,ip);
    estimate = exp(logEstimate-max(logEstimate));
    estimate = estimate/trapz(y_grid, estimate);

    errors = squeeze(errorDraws(:,ip,:));
    logDraws = logReference + shockDeltaX(ip)*(est_all(:,ip)-errors);
    logDraws = logDraws-max(logDraws, [], 1);
    densityDraws = exp(logDraws);
    areas = trapz(y_grid, densityDraws, 1);
    if any(~isfinite(areas)) || any(areas <= 0)
        error('A simulated implied density could not be normalized.');
    end
    densityDraws = densityDraws./areas;
    bounds = quantile(densityDraws, [alpha, 1-alpha], 2);

    DensityCI(ip).estimate = estimate;
    DensityCI(ip).lower = bounds(:,1);
    DensityCI(ip).upper = bounds(:,2);
end
end
