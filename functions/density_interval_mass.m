function p = density_interval_mass(x, f, lowerBound, upperBound)
%DENSITY_INTERVAL_MASS Probability over an interval under the grid CDF.

x = x(:);
f = max(f(:),0);
area = trapz(x,f);
if numel(x) ~= numel(f) || any(diff(x) <= 0) || ~(area > 0)
    error('x and f must define a nondegenerate density on an increasing grid.');
end
f = f / area;
cdf = cumtrapz(x,f);
cdf = cdf / cdf(end);

lowerProbability = local_cdf(lowerBound,x,cdf);
upperProbability = local_cdf(upperBound,x,cdf);
p = max(0,min(1,upperProbability-lowerProbability));
end

function value = local_cdf(point,x,cdf)
if point <= x(1)
    value = 0;
elseif point >= x(end)
    value = 1;
else
    value = interp1(x,cdf,point,'linear');
end
end
