function p = normalize_density_trapz(f, grid)
%NORMALIZE_DENSITY_TRAPZ Normalize a density by trapezoidal integration.
    f = f(:);
    grid = grid(:);

    f(~isfinite(f)) = 0;
    f(f < 0) = 0;

    area = trapz(grid, f);
    if area <= 0 || ~isfinite(area)
        span = grid(end) - grid(1);
        if span <= 0
            p = ones(size(f)) / numel(f);
        else
            p = ones(size(f)) / span;
        end
    else
        p = f / area;
    end
end
