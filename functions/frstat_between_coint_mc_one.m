function out = frstat_between_coint_mc_one(U_pre, xvar, version, varargin)
%FRSTAT_BETWEEN_COINT_MC_ONE Between-cointegration Monte Carlo test.
%
% version = 1: no covariance denominator
% version = 2: covariance denominator along the dominant residual-covariance
%              direction.

p = inputParser;
addRequired(p, 'U_pre', @(x) isnumeric(x) && ~isempty(x));
addRequired(p, 'xvar', @(x) isnumeric(x) && ~isempty(x));
addRequired(p, 'version', @(x) isnumeric(x) && isscalar(x) && any(x == [1 2]));
addParameter(p, 'Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'Bmc', 5000, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'McGrid', 499, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'TestDemean', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Seed', 12345, @(x) isnumeric(x) && isscalar(x));
parse(p, U_pre, xvar, version, varargin{:});

alpha = p.Results.Alpha;
Bmc = p.Results.Bmc;
mcgrid = p.Results.McGrid;
test_demean = p.Results.TestDemean;
seed = p.Results.Seed;

oldRng = rng;
rngCleanup = onCleanup(@() rng(oldRng));
rng(seed, 'twister');

TTT = size(U_pre, 2);
bandwidth = round(TTT^(1/4));

Xtest = double(xvar);
zz = double(U_pre.');

if test_demean
    Xtest = Xtest - mean(Xtest, 1);
    zz = zz - mean(zz, 1);
end

TTT = size(zz, 1);
KX = size(Xtest, 2);
KU = size(zz, 2);

if size(Xtest,1) ~= TTT
    error('Xtest and residual series zz have different time lengths.');
end

zz2 = cumsum(zz, 1);
vx1 = (zz2.' * zz2) / (TTT^2);
vx1 = (vx1 + vx1.') / 2;

if version == 1
    teststat = max(real(eig(vx1)));
    ev0 = [];
    eval_vxcov = [];
else
    vxcov = (zz.' * zz) / TTT;
    vxcov = (vxcov + vxcov.') / 2;
    [Vcov, Dcov] = eig(vxcov);
    [~, idx] = sort(real(diag(Dcov)), 'descend');
    ev0 = real(Vcov(:,idx(1)));

    eval_vx1 = ev0.' * vx1 * ev0;
    eval_vxcov = ev0.' * vxcov * ev0;
    if ~(isfinite(eval_vxcov) && eval_vxcov > 1e-10)
        error('Residual covariance denominator is too small or non-finite.');
    end
    teststat = eval_vx1 / eval_vxcov;
end

XDtest = Xtest(2:TTT,:) - Xtest(1:(TTT-1),:);
UDtest = zz(2:TTT,:);
Ztest = [XDtest, UDtest];

OmegaZZhat = lr_var(Ztest,Ztest,2,bandwidth).omega;
OmegaZZhat = (OmegaZZhat + OmegaZZhat.') / 2;

GammaXUplus = lr_var(UDtest,XDtest,1,bandwidth).omega;
if ~isequal(size(GammaXUplus), [KU, KX])
    if isequal(size(GammaXUplus), [KX, KU])
        GammaXUplus = GammaXUplus.';
    else
        error('Unexpected GammaXUplus dimensions.');
    end
end

sqrtOmegaZZ = frstat_sym_psd_sqrt(OmegaZZhat);
dt_mc = 1 / mcgrid;
grid_mc = linspace(0, 1, mcgrid + 1);

T0mc = NaN(Bmc, 1);
for bb = 1:Bmc
    dW = sqrtOmegaZZ * randn(KX + KU, mcgrid) * sqrt(dt_mc);

    Wpath = zeros(KX + KU, mcgrid + 1);
    Wpath(:,2:end) = cumsum(dW, 2);

    Wx = Wpath(1:KX,:);
    Wu = Wpath((KX+1):(KX+KU),:);
    dWu = dW((KX+1):(KX+KU),:);

    if test_demean
        Wx_left0 = Wx(:,1:mcgrid);
        Wx_bar = sum(Wx_left0, 2) * dt_mc;
        Wx_use = Wx - Wx_bar * ones(1, mcgrid + 1);
        Wu_use = Wu - Wu(:,mcgrid + 1) * grid_mc;
    else
        Wx_use = Wx;
        Wu_use = Wu;
    end

    Wx_left = Wx_use(:,1:mcgrid);
    Jhat = Wx_left * Wx_left.' * dt_mc;
    Jhat = (Jhat + Jhat.') / 2;
    Jinv = (Jhat + eye(KX) * 1e-10) \ eye(KX);

    Qhat = zeros(KX, mcgrid + 1);
    Qhat(:,2:end) = cumsum(Wx_left * dt_mc, 2);

    Ahat = dWu * Wx_left.' + GammaXUplus;
    Rhat = Wu_use - Ahat * Jinv * Qhat;

    Rleft = Rhat(:,1:mcgrid);
    KhatR = Rleft * Rleft.' * dt_mc;
    KhatR = (KhatR + KhatR.') / 2;

    if version == 1
        numerR = max(real(eig(KhatR)));
        if isfinite(numerR)
            T0mc(bb) = numerR;
        end
    else
        numerR = ev0.' * KhatR * ev0;
        denomR = eval_vxcov;
        if isfinite(numerR) && isfinite(denomR) && denomR > 1e-10
            T0mc(bb) = numerR / denomR;
        end
    end
end

T0mc = T0mc(isfinite(T0mc));
if numel(T0mc) ~= Bmc
    error('frstat_between_coint_mc_one:InvalidMonteCarloDraws', ...
        'Only %d of %d requested Monte Carlo draws were finite.', ...
        numel(T0mc), Bmc);
end
qu = frstat_quantile_type7(T0mc, 1 - alpha);
pval = (1 + sum(T0mc >= teststat)) / (numel(T0mc) + 1);

alphas = [0.10 0.05 0.01];
critvals = array2table(NaN(1, numel(alphas)), ...
    'VariableNames', {'alpha_0_10','alpha_0_05','alpha_0_01'});
for ii = 1:numel(alphas)
    critvals{1,ii} = frstat_quantile_type7(T0mc, 1 - alphas(ii));
end

if pval > alpha
    decision = "Between-cointegrated: fail to reject H0";
else
    decision = "Reject between-cointegration";
end

out = struct('version', version, 'alpha', alpha, 'Bmc', Bmc, ...
    'mcgrid', mcgrid, 'test_demean', test_demean, 'seed', seed, ...
    'teststat', teststat, 'qu', qu, 'pval', pval, ...
    'nvalid', numel(T0mc), 'T0mc', T0mc, 'critvals', critvals, ...
    'decision', decision, 'ev0', ev0, 'eval_vxcov', eval_vxcov, ...
    'OmegaZZhat', OmegaZZhat, 'GammaXUplus', GammaXUplus);
end
