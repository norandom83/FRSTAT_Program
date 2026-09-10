function out = lr_var(u, v, side, bandwidth)
%LR_VAR One- or two-sided long-run covariance using the Parzen kernel.

u = double(u);
v = double(v);
if size(u,1) ~= size(v,1)
    error('u and v must have the same number of observations.');
end

T = size(u,1);
nLag = floor(bandwidth);
lagRatio = (1:nLag)'/bandwidth;
weights = (1-6*lagRatio.^2+6*lagRatio.^3).*(lagRatio <= 0.5) + ...
    2*(1-lagRatio).^3.*(lagRatio > 0.5);

omega = u'*v/T;
if side == 2
    for lag = 1:nLag
        omega = omega + weights(lag)*( ...
            u(1:T-lag,:)'*v(1+lag:T,:) + ...
            u(1+lag:T,:)'*v(1:T-lag,:))/T;
    end
elseif side == 1
    for lag = 1:nLag
        omega = omega + weights(lag)*u(1+lag:T,:)'*v(1:T-lag,:)/T;
    end
else
    error('side must be 1 or 2.');
end

out = struct('omega',omega,'bandw',bandwidth);
end
