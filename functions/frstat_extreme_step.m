function out = frstat_extreme_step(y_grid, f_ref0, est, deltaX, locShift, S, forcing)
%FRSTAT_EXTREME_STEP Cold, warm, and total extreme-state probability changes.

y_grid = y_grid(:);
[f_ref, f_full] = clr_push_density_from_est(y_grid, f_ref0, est, deltaX);
f_loc = shift_density_on_grid(y_grid, f_ref, locShift);
rLow = density_quantile_from_grid(y_grid, f_ref, S.q_low);
rHigh = density_quantile_from_grid(y_grid, f_ref, S.q_high);
coldRef = density_interval_mass(y_grid,f_ref,-Inf,rLow);
warmRef = density_interval_mass(y_grid,f_ref,rHigh,Inf);
dCold = density_interval_mass(y_grid,f_full,-Inf,rLow)-coldRef;
dWarm = density_interval_mass(y_grid,f_full,rHigh,Inf)-warmRef;
dColdLoc = density_interval_mass(y_grid,f_loc,-Inf,rLow)-coldRef;
dWarmLoc = density_interval_mass(y_grid,f_loc,rHigh,Inf)-warmRef;

out = struct('forcing',string(forcing),'q_low',S.q_low,'q_high',S.q_high, ...
    'r_low',rLow,'r_high',rHigh,'cold_ref',coldRef,'warm_ref',warmRef, ...
    'total_ref',coldRef+warmRef,'d_cold',dCold,'d_warm',dWarm, ...
    'd_total',dCold+dWarm,'d_cold_loc',dColdLoc,'d_warm_loc',dWarmLoc, ...
    'd_total_loc',dColdLoc+dWarmLoc,'d_cold_dist',dCold-dColdLoc, ...
    'd_warm_dist',dWarm-dWarmLoc, ...
    'd_total_dist',(dCold+dWarm)-(dColdLoc+dWarmLoc));
end
