function Ahalf = frstat_sym_psd_sqrt(S, tol)
%FRSTAT_SYM_PSD_SQRT Symmetric positive semidefinite square root.

if nargin < 2
    tol = 1e-10;
end

S = (S + S.') / 2;
[V, D] = eig(S);
vals = real(diag(D));
scalev = max(max(abs(vals)), 1);
vals = max(vals, tol * scalev);

Ahalf = real(V) * diag(sqrt(vals)) * real(V).';
end
