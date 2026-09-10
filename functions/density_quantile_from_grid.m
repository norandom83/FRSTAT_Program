function qx = density_quantile_from_grid(x, f, q)
% Quantile of a density f on grid x using CDF interpolation (robust to flat CDF segments).

    x = x(:); f = f(:);
    f = max(f,0);

    area = trapz(x,f);
    if ~(area > 0)
        qx = NaN;
        return;
    end
    f = f / area;

    cdf = cumtrapz(x,f);
    cdf = max(min(cdf,1),0);

    % --- FIX: interp1 requires unique sample points ---
    % Keep only the first occurrence of each unique CDF value (stable).
    [cdf_u, ia] = unique(cdf, 'stable');
    x_u = x(ia);

    % Guard: if CDF is degenerate (all same), quantile undefined
    if numel(cdf_u) < 2
        qx = NaN;
        return;
    end

    % Clamp q within [min,max] to avoid extrap explosions
    q = min(max(q, cdf_u(1)), cdf_u(end));

    qx = interp1(cdf_u, x_u, q, 'linear');
end
