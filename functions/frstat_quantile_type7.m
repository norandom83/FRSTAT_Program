function q = frstat_quantile_type7(x, p)
%FRSTAT_QUANTILE_TYPE7 R quantile(..., type = 7) for deterministic parity.

x = sort(x(:));
x = x(isfinite(x));
n = numel(x);

if n == 0
    q = NaN(size(p));
    return;
end

origSize = size(p);
p = p(:).';

h = (n - 1) .* p + 1;
lo = floor(h);
hi = ceil(h);
lo = max(1, min(n, lo));
hi = max(1, min(n, hi));
gamma = h - lo;

q = (1 - gamma) .* x(lo).' + gamma .* x(hi).';
q = reshape(q, origSize);
end
