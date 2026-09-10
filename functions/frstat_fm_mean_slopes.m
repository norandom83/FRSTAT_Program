function out = frstat_fm_mean_slopes(y, x)
%FRSTAT_FM_MEAN_SLOPES Fully modified slopes for a scalar response.

y = double(y(:));
x = double(x);
if size(x,1) ~= numel(y) || any(~isfinite(y)) || any(~isfinite(x(:)))
    error('y and x must be finite and have the same number of observations.');
end
if size(x,1) <= size(x,2) + 2
    error('The scalar benchmark has too few observations for its predictors.');
end

yCentered = y - mean(y);
xCentered = x - mean(x,1);
preliminarySlope = (xCentered' * xCentered) \ (xCentered' * yCentered);
preliminaryResidual = yCentered - xCentered * preliminarySlope;

deltaX = diff(xCentered,1,1);
residual = preliminaryResidual(2:end,:);
kerband = round(size(deltaX,1)^(1/4));

omegaXX = lr_var(deltaX,deltaX,2,kerband).omega;
omegaXX = (omegaXX+omegaXX')/2;
if rcond(omegaXX)<1e-12
    error('The scalar-benchmark forcing covariance is numerically singular.');
end
omegaUX = lr_var(residual,deltaX,2,kerband).omega;
omegaUXplus = lr_var(residual,deltaX,1,kerband).omega;
omegaXXplus = lr_var(deltaX,deltaX,1,kerband).omega;
upsilon = omegaUXplus - omegaUX / omegaXX * omegaXXplus;

xEst = xCentered(2:end,:);
xEst = xEst - mean(xEst,1);
yEst = yCentered(2:end)';
zEst = yEst - omegaUX / omegaXX * deltaX';
nEst = size(xEst,1);
slope = ((xEst' * xEst) \ (xEst' * zEst' - nEst * upsilon'))';

out = struct();
out.slope = slope;
out.intercept = mean(y) - slope * mean(x,1)';
out.preliminarySlope = preliminarySlope';
out.preliminaryResidual = preliminaryResidual;
out.nEst = nEst;
out.bandwidth = kerband;
out.upsilon = upsilon;
end
