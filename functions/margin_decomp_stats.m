function mstats = margin_decomp_stats(y_grid, f_ref, f_new, q_high)
% Compute median / IQR / upper-tail mass differences from two densities.
% Robust to flat CDF regions.

    y_grid = y_grid(:);

    f_ref = max(f_ref(:), 0);
    f_new = max(f_new(:), 0);

    Aref = trapz(y_grid, f_ref);
    Anew = trapz(y_grid, f_new);

    if Aref <= 0 || Anew <= 0
        error('Density has non-positive integral.');
    end

    f_ref = f_ref / Aref;
    f_new = f_new / Anew;

    cdf_ref = cumtrapz(y_grid, f_ref);
    cdf_new = cumtrapz(y_grid, f_new);

    cdf_ref = cdf_ref / max(cdf_ref(end), eps);
    cdf_new = cdf_new / max(cdf_new(end), eps);

    cdf_ref = cummax(cdf_ref);
    cdf_new = cummax(cdf_new);

    q25_ref = invcdf_monotone(cdf_ref, y_grid, 0.25);
    q50_ref = invcdf_monotone(cdf_ref, y_grid, 0.50);
    q75_ref = invcdf_monotone(cdf_ref, y_grid, 0.75);

    q25_new = invcdf_monotone(cdf_new, y_grid, 0.25);
    q50_new = invcdf_monotone(cdf_new, y_grid, 0.50);
    q75_new = invcdf_monotone(cdf_new, y_grid, 0.75);

    r_high = invcdf_monotone(cdf_ref, y_grid, q_high);

    tail_ref = density_interval_mass(y_grid, f_ref, r_high, Inf);
    tail_new = density_interval_mass(y_grid, f_new, r_high, Inf);

    mstats.r_high     = r_high;

    mstats.median_ref = q50_ref;
    mstats.median_new = q50_new;
    mstats.d_median   = q50_new - q50_ref;

    mstats.iqr_ref    = q75_ref - q25_ref;
    mstats.iqr_new    = q75_new - q25_new;
    mstats.d_iqr      = mstats.iqr_new - mstats.iqr_ref;

    mstats.tail_ref   = tail_ref;
    mstats.tail_new   = tail_new;
    mstats.d_tail     = tail_new - tail_ref;
end

function yq = invcdf_monotone(cdf, y_grid, q)
    q = min(max(q, 0), 1);

    [cdf_u, ia] = unique(cdf, 'last');
    y_u = y_grid(ia);

    if cdf_u(1) > 0
        cdf_u = [0; cdf_u];
        y_u   = [y_grid(1); y_u];
    end
    if cdf_u(end) < 1
        cdf_u = [cdf_u; 1];
        y_u   = [y_u; y_grid(end)];
    end

    yq = interp1(cdf_u, y_u, q, 'linear');
end
