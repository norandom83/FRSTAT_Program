function [f_ref, f_new] = clr_push_density_from_est(y_grid, yvar0, est, deltaX)
% Build reference density f_ref on y_grid and apply CLR push using est and deltaX.
% est is a CLR response function (length(y_grid) x 1), deltaX is scalar shock size.

    y_grid = y_grid(:);
    est    = est(:);

    % Reference density: if matrix, take mean across columns; if vector, use directly
    if ismatrix(yvar0) && size(yvar0,2) > 1
        f_ref = mean(yvar0, 2, 'omitnan');
    else
        f_ref = yvar0(:);
    end

    % Ensure nonnegativity and normalization
    f_ref = max(f_ref, 0);
    area0 = trapz(y_grid, f_ref);
    if area0 <= 0
        error('Reference density has non-positive integral. Check yvar0/y_grid.');
    end
    f_ref = f_ref / area0;

    % CLR of reference using the same trapezoidal integral as the estimator.
    eps0 = 1e-12;
    logf0 = log(max(f_ref, eps0));
    mlogf0 = trapz(y_grid, logf0) / (y_grid(end)-y_grid(1));
    clr0 = logf0 - mlogf0;

    % Apply shock in CLR space
    clr1 = clr0 + deltaX * est;

    % Invert CLR -> density and renormalize
    f_new = exp(clr1);
    f_new = max(f_new, 0);
    area1 = trapz(y_grid, f_new);
    f_new = f_new / area1;
end
