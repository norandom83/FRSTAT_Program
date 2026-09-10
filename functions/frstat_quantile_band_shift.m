function delta_loc = frstat_quantile_band_shift(y_grid, f_ref, f_full, u_band, u_n)
% Compute central-quantile average displacement:
% delta_loc = average_u [ Q_full(u) - Q_ref(u) ], u in [uL, uH]

    if nargin < 5 || isempty(u_n)
        u_n = 201;
    end

    y_grid = y_grid(:);
    f_ref  = max(f_ref(:),  0);
    f_full = max(f_full(:), 0);

    f_ref  = f_ref  / trapz(y_grid, f_ref);
    f_full = f_full / trapz(y_grid, f_full);

    uL = u_band(1);
    uH = u_band(2);

    u_grid = linspace(uL, uH, u_n)';

    q_ref  = NaN(u_n,1);
    q_full = NaN(u_n,1);

    for i = 1:u_n
        q_ref(i)  = density_quantile_from_grid(y_grid, f_ref,  u_grid(i));
        q_full(i) = density_quantile_from_grid(y_grid, f_full, u_grid(i));
    end

    delta_loc = trapz(u_grid, q_full - q_ref) / (uH - uL);
end