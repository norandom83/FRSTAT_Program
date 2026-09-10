function out = frstat_local_loss_ci(y_grid, f_ref0, est_all, betaDraws, deltaX, S, localProfiles)
%FRSTAT_LOCAL_LOSS_CI Pointwise bands for local probability-mass responses.

y_grid = y_grid(:);
f_ref = build_reference_density_seed(y_grid, f_ref0);
h = S.loss_h;
fineStep = min(h/20, mean(diff(y_grid))/5);
sFine = (y_grid(1):fineStep:y_grid(end))';
if sFine(end) < y_grid(end)
    sFine = [sFine; y_grid(end)];
end
rGrid = localProfiles(1).r_grid(:);
weights = local_window_weights(sFine, rGrid, h);
fRefFine = normalize_density_trapz( ...
    interp1(y_grid, f_ref, sFine, 'linear', 0), sFine);

nForce = size(est_all, 2);
nRep = size(betaDraws, 3);
bands = repmat(struct('lower', [], 'upper', []), nForce, 1);
for ip = 1:nForce
    profiles = NaN(numel(rGrid), nRep);
    for b = 1:nRep
        % betaDraws simulate estimation error (estimate minus truth).
        betaDraw = est_all(:,ip) - betaDraws(:,ip,b);
        f = exp(log(max(f_ref, realmin)) + deltaX(ip) * betaDraw);
        f = normalize_density_trapz(f, y_grid);
        fFine = normalize_density_trapz( ...
            interp1(y_grid, f, sFine, 'linear', 0), sFine);
        profiles(:,b) = weights * (fFine - fRefFine);
    end
    bands(ip).lower = quantile(profiles, S.ciAlpha, 2);
    bands(ip).upper = quantile(profiles, 1-S.ciAlpha, 2);
end
out = struct('r_grid', rGrid, 'bands', bands, ...
    'nRep', nRep, 'alpha', S.ciAlpha);
end

function W = local_window_weights(grid, centers, width)
W = zeros(numel(centers), numel(grid));
for i = 1:numel(centers)
    idx = find(grid >= centers(i)-width/2 & grid <= centers(i)+width/2);
    if numel(idx) < 2
        continue;
    end
    x = grid(idx);
    w = zeros(numel(idx),1);
    w(1) = (x(2)-x(1))/2;
    w(end) = (x(end)-x(end-1))/2;
    if numel(idx) > 2
        w(2:end-1) = (x(3:end)-x(1:end-2))/2;
    end
    W(i,idx) = w;
end
end
