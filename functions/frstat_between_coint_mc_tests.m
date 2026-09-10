function out = frstat_between_coint_mc_tests(U_pre, xvar, varargin)
%FRSTAT_BETWEEN_COINT_MC_TESTS Run the two between-cointegration diagnostics.

p = inputParser;
addRequired(p, 'U_pre', @(x) isnumeric(x) && ~isempty(x));
addRequired(p, 'xvar', @(x) isnumeric(x) && ~isempty(x));
addParameter(p, 'Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'Bmc', 5000, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'McGrid', 499, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'TestDemean', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Seed', 12345, @(x) isnumeric(x) && isscalar(x));
parse(p, U_pre, xvar, varargin{:});

args = {'Alpha', p.Results.Alpha, 'Bmc', p.Results.Bmc, ...
    'McGrid', p.Results.McGrid, 'TestDemean', p.Results.TestDemean, ...
    'Seed', p.Results.Seed};

out = struct();
out.version1 = frstat_between_coint_mc_one(U_pre, xvar, 1, args{:});
out.version2 = frstat_between_coint_mc_one(U_pre, xvar, 2, args{:});
out.alpha = p.Results.Alpha;
out.Bmc = p.Results.Bmc;
out.mcgrid = p.Results.McGrid;
out.test_demean = p.Results.TestDemean;
out.seed = p.Results.Seed;
end
