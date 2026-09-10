function out = frstat_local_loss_profile_step(y_grid, f_ref0, est, deltaX, loc_shift, h, step, forcing)

    y_grid = y_grid(:);
    f_ref0 = f_ref0(:);
    est    = est(:);

    % ---- full counterfactual density from CLR push ----
    eps0 = 1e-12;
    g_full = log(max(f_ref0, eps0)) + deltaX * est;
    f_full = exp(g_full);
    f_full = f_full / trapz(y_grid, f_full);

    % ---- location-only counterfactual density ----
    % f_loc(s) = f0(s - loc_shift)
    f_loc = interp1(y_grid, f_ref0, y_grid - loc_shift, 'linear', 0);
    area_loc = trapz(y_grid, f_loc);
    if area_loc <= 0
        error('Location-only density has non-positive integral.');
    end
    f_loc = f_loc / area_loc;

    % =====================================================
    % Use a finer grid so that h=0.10 bins are meaningfully integrated
    % =====================================================
    dy = mean(diff(y_grid));
    fineStep = min(h/20, dy/5);   % sufficiently finer than both h and original grid
    s_fine = (y_grid(1):fineStep:y_grid(end))';
    if s_fine(end) < y_grid(end)
        s_fine = [s_fine; y_grid(end)];
    end

    f_ref0_f = interp1(y_grid, f_ref0, s_fine, 'linear', 0);
    f_full_f = interp1(y_grid, f_full, s_fine, 'linear', 0);
    f_loc_f  = interp1(y_grid, f_loc,  s_fine, 'linear', 0);

    % re-normalize on fine grid for numerical stability
    f_ref0_f = f_ref0_f / trapz(s_fine, f_ref0_f);
    f_full_f = f_full_f / trapz(s_fine, f_full_f);
    f_loc_f  = f_loc_f  / trapz(s_fine, f_loc_f);

   % =====================================================
% Overlapping rolling windows of width h across the support
% =====================================================
% step = 0.02;   % overlap step; smaller than h
r_grid = (y_grid(1)+h/2 : step : y_grid(end)-h/2)';   % window centers

nR = numel(r_grid);
dLoss_full = zeros(nR,1);
dLoss_loc  = zeros(nR,1);

for i = 1:nR
    lo = r_grid(i) - h/2;
    hi = r_grid(i) + h/2;

    idx = (s_fine >= lo & s_fine <= hi);

    dLoss_full(i) = trapz(s_fine(idx), f_full_f(idx) - f_ref0_f(idx));
    dLoss_loc(i)  = trapz(s_fine(idx), f_loc_f(idx)  - f_ref0_f(idx));
end

    out.forcing    = char(forcing);
    out.r_grid     = r_grid;
    out.dLoss_full = dLoss_full;
    out.dLoss_loc  = dLoss_loc;
    out.dLoss_dist = dLoss_full - dLoss_loc;
end