function f_shift = shift_density_on_grid(y_grid, f_ref, dmu)
% Mean-only counterfactual: shift density by dmu using interpolation on the same grid.

    y_grid = y_grid(:);
    f_ref = f_ref(:);
    f_ref = max(f_ref,0);
    f_ref = f_ref / trapz(y_grid,f_ref);

    % f_shift(y) = f_ref(y - dmu)
    f_shift = interp1(y_grid, f_ref, y_grid - dmu, 'linear', 0);
    f_shift = max(f_shift,0);
    f_shift = f_shift / trapz(y_grid,f_shift);
end
