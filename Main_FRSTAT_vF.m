%% Main FRSTAT vF Program
%
%                                                           
%                                                           

clear; clc; close all;

codeDir = fileparts(mfilename('fullpath'));
addpath(fullfile(codeDir, 'functions'));
dataFile = fullfile(codeDir, 'FRSTAT_Data.xlsx');

seed_response = 12345;
seed_test = 12345;

%% Data
ForcingInput = readmatrix(dataFile, 'Sheet', 'Forcing_Components', ...
    'Range', 'A2:L176');
forcingYears = ForcingInput(:,1);
xvar = [ForcingInput(:,2), sum(ForcingInput(:,3:end),2)];
Fnames = ["F1", "F2"];
forcingCentered = xvar-mean(xvar,1);
F1F2Correlation = (forcingCentered(:,1)'*forcingCentered(:,2))/ ...
    sqrt(sum(forcingCentered(:,1).^2)*sum(forcingCentered(:,2).^2));

BaselineInput = readmatrix(dataFile, 'Sheet', 'Baseline_Density', 'Range', 'A2');
HemisphericInput = readmatrix(dataFile, 'Sheet', 'Hemispheric_Density', 'Range', 'A2');
FixedGridInput = readmatrix(dataFile, 'Sheet', 'Fixed_Grid_Density', 'Range', 'A2');
SampleInfo = readtable(dataFile, 'Sheet', 'Sample_Info', ...
    'Range', 'A1:B5', 'VariableNamingRule', 'preserve');
CoverageInput = readmatrix(dataFile, 'Sheet', 'Coverage_Annual', ...
    'Range', 'A2:D176');

y_grid = BaselineInput(:,1);
yvar0 = BaselineInput(:,2:end);
y_grid_hemi = HemisphericInput(:,1);
yvar0_hemi = HemisphericInput(:,2:end);
y_grid_fixed = FixedGridInput(:,1);
yvar0_fixed = FixedGridInput(:,2:end);

sampleItems = string(SampleInfo.Item);
sampleStartYear = SampleInfo.Value(sampleItems == "StartYear");
sampleEndYear = SampleInfo.Value(sampleItems == "EndYear");
totalGridCells = SampleInfo.Value(sampleItems == "TotalGridCellCount");
fixedGridCells = SampleInfo.Value(sampleItems == "FixedGridCellCount");
expectedYears = (sampleStartYear:sampleEndYear)';
coverageYears = CoverageInput(:,1);
observedCells = CoverageInput(:,2);
observedCellMonths = CoverageInput(:,3);
fixedGridCellMonths = CoverageInput(:,4);
fixedGridShare = fixedGridCells/totalGridCells;
medianObservedCellMonths = median(observedCellMonths);
medianFixedGridCellMonths = median(fixedGridCellMonths);
if ~isequal(forcingYears,expectedYears) || ~isequal(coverageYears,expectedYears)
    error('Forcing, density, and coverage years are not aligned.');
end

%% CLR transformation and basis representation
nt = numel(y_grid);
Yraw = NaN(size(yvar0));
for i = 1:size(yvar0,2)
    logDensity = log(yvar0(:,i));
    clrMean = trapz(y_grid, logDensity) / (y_grid(end) - y_grid(1));
    Yraw(:,i) = logDensity - clrMean;
end

t = (0:(nt-1))'/(nt-1);
inner = @(f,g) trapz(t, f.*g);
lbnumber2 = 20;
LBF = NaN(nt, lbnumber2);

for i = 1:(lbnumber2/2)
    s = sqrt(2)*sin(2*pi*i*t);
    c = sqrt(2)*cos(2*pi*i*t);
    LBF(:,2*i-1) = s / sqrt(inner(s,s));
    LBF(:,2*i) = c / sqrt(inner(c,c));
end

LBF = [ones(nt,1), LBF];
LBF(:,1) = LBF(:,1) / sqrt(inner(LBF(:,1),LBF(:,1)));
for i = 2:size(LBF,2)
    for j = 1:size(LBF,2)
        if j ~= i
            LBF(:,i) = LBF(:,i) - ...
                inner(LBF(:,i),LBF(:,j))/inner(LBF(:,j),LBF(:,j))*LBF(:,j);
        end
    end
end
for i = 1:size(LBF,2)
    LBF(:,i) = LBF(:,i) / sqrt(inner(LBF(:,i),LBF(:,i)));
end

trapWeights = ones(nt,1)/(nt-1);
trapWeights([1,end]) = trapWeights([1,end])/2;
YY = LBF' * (Yraw .* trapWeights);
YY = YY - mean(YY,2);
T = size(YY,2);

%% Fully modified estimation
XX = xvar - mean(xvar,1);
beta_pre = ((XX' * XX) \ (XX' * YY'))';
U_pre = YY - beta_pre * XX';

XD = diff(XX);
UD = U_pre(:,2:end)';
kerband = round(size(XD,1)^(1/4));

OmegaXX2 = lr_var(XD,XD,2,kerband).omega;
OmegaXX2 = (OmegaXX2 + OmegaXX2')/2;
if rcond(OmegaXX2) < 1e-12
    error('The forcing long-run covariance is numerically singular.');
end
OmegaXU2 = lr_var(UD,XD,2,kerband).omega;
OmegaXU1 = lr_var(UD,XD,1,kerband).omega;
OmegaXX1 = lr_var(XD,XD,1,kerband).omega;
Upsilon = OmegaXU1 - OmegaXU2/OmegaXX2*OmegaXX1;

XX2 = XX(2:end,:);
XX2 = XX2 - mean(XX2,1);
Z1 = YY(:,2:end) - OmegaXU2/OmegaXX2*XD';
nEst = size(XX2,1);

% The one-sided correction is N*Upsilon in cross-product form.
beta = ((XX2' * XX2) \ (XX2' * Z1' - nEst*Upsilon'))';

OmegaUU2 = lr_var(UD,UD,2,kerband).omega;
OmegaUN = OmegaUU2 - OmegaXU2/OmegaXX2*OmegaXU2';
OmegaUN = (OmegaUN + OmegaUN')/2;
omegaUNTolerance = 1e-10*max(1,norm(OmegaUN,2));
if min(eig(OmegaUN)) < -omegaUNTolerance
    error('The estimated conditional response covariance is materially non-PSD.');
end

KK = round(sqrt(T/kerband));
[VecUN,DUN] = eig(OmegaUN,'vector');
[evals,idxUN] = sort(DUN,'descend');
evecs = VecUN(:,idxUN);
evals = evals(1:KK);
evecs = evecs(:,1:KK);

[VecX,DX] = eig(OmegaXX2,'vector');
[evalsx,idxX] = sort(DX,'descend');
evecsx = VecX(:,idxX);
if min(evalsx) <= 0
    error('The forcing long-run covariance must be positive definite.');
end
evecs = LBF*evecs;

%% Functional responses and margins
S = frstat_default_settings();
nForce = size(xvar,2);
ShockDeltaX = std(xvar,0,1)';              % Positive one-s.d. shocks
restrictionReplications = 50000;

est_all = NaN(numel(y_grid),nForce);
CIband = repmat(struct('estimate',[],'lquan',[],'uquan',[]),nForce,1);
f_ref0_loss = build_reference_density_seed(y_grid,yvar0);

% Joint draws are shared across response contrasts and nonlinear margins.
RestrictionDraws = frstat_simulate_beta_draws(nForce,KK,S.bmgrid,XX2, ...
    evecs,evals,evecsx,evalsx,nEst,restrictionReplications,seed_response);
if size(RestrictionDraws,3) ~= restrictionReplications || ...
        any(~isfinite(RestrictionDraws),'all')
    error('The joint coefficient simulation did not return all requested draws.');
end
BetaDraws = RestrictionDraws(:,:,1:S.nRep);

for ip = 1:nForce
    est = (beta(:,ip)'*LBF')';
    [lquan,uquan,localProjection] = compute_pointwise_ci( ...
        squeeze(BetaDraws(:,ip,:)),y_grid,S.stepsize);

    est_all(:,ip) = est;
    CIband(ip).estimate = localProjection*est;
    CIband(ip).lquan = lquan;
    CIband(ip).uquan = uquan;

    deltaX = ShockDeltaX(ip);
    [tmpMarginStats,tmpMarginDecomp] = frstat_margin_step( ...
        y_grid,yvar0,est,deltaX,S,Fnames(ip));
    tmpLocalLossCurve = frstat_local_loss_profile_step( ...
        y_grid,f_ref0_loss,est,deltaX,tmpMarginDecomp.loc_shift, ...
        S.loss_h,S.loss_step,Fnames(ip));
    tmpExtremeStats = frstat_extreme_step( ...
        y_grid,f_ref0_loss,est,deltaX,tmpMarginDecomp.loc_shift,S,Fnames(ip));

    if ip == 1
        MarginStats = repmat(tmpMarginStats,nForce,1);
        MarginDecomp = repmat(tmpMarginDecomp,nForce,1);
        LocalLossCurve = repmat(tmpLocalLossCurve,nForce,1);
        ExtremeStats = repmat(tmpExtremeStats,nForce,1);
    end
    MarginStats(ip) = tmpMarginStats;
    MarginDecomp(ip) = tmpMarginDecomp;
    LocalLossCurve(ip) = tmpLocalLossCurve;
    ExtremeStats(ip) = tmpExtremeStats;
end

LocalLossCI = frstat_local_loss_ci(y_grid,f_ref0_loss,est_all,RestrictionDraws, ...
    ShockDeltaX,S,LocalLossCurve);
DensityCI = frstat_density_ci(y_grid,f_ref0_loss,est_all,RestrictionDraws, ...
    ShockDeltaX,S.ciAlpha);
CommonResponseTest = frstat_common_response_test( ...
    est_all,y_grid,evecs,evals,evecsx,evalsx,XX2,nEst, ...
    S.bmgrid,restrictionReplications,seed_response,RestrictionDraws);
clear RestrictionDraws

CommonResponseSensitivity = frstat_common_response_sensitivity( ...
    y_grid,yvar0,xvar,S,restrictionReplications,seed_response, ...
    CommonResponseTest);

%% Cointegration diagnostics and coverage checks
CVInput = readmatrix(dataFile, 'Sheet', 'Breitung_CV', 'Range', 'A2');
CVx = CVInput(:,2:end);
WithinCointTest = frstat_within_coint_test(xvar,CVx,0.05);

BetweenCointMC = frstat_between_coint_mc_tests( ...
    U_pre,xvar,'Alpha',0.05,'Bmc',S.testReplications, ...
    'McGrid',S.testMcGrid,'TestDemean',true,'Seed',seed_test);

HemisphericCheck = frstat_reestimate_check( ...
    y_grid_hemi,yvar0_hemi,xvar,S,lbnumber2,seed_test, ...
    "hemispheric equal weight");
FixedGridCheck = frstat_reestimate_check( ...
    y_grid_fixed,yvar0_fixed,xvar,S,lbnumber2,seed_test,"fixed grid");

%% Tables
ForcingAggregation = table( ...
    ["F1: CO2 forcing"; "F2: Non-CO2 anthropogenic forcing"], ...
    ["CO2"; ...
     "CH4, N2O, aerosol-radiation interactions, aerosol-cloud interactions, O3, contrails, land-use change, BC on snow, H2O_strat, halogenated species"], ...
    ["The dominant anthropogenic greenhouse-gas component in the IPCC effective radiative forcing accounting. This portfolio provides the main low-frequency CO2 forcing benchmark."; ...
     "Anthropogenic forcing outside the CO2 block. This portfolio collects non-CO2 greenhouse gases, reactive-chemistry, aerosol, cloud-adjustment, surface-albedo, snow-albedo, and aviation-related forcing channels."], ...
    'VariableNames',{'Group','ForcingComponents','StatisticalRole'});

PaperMarginTable = table(Fnames(:),[MarginStats.dmu_full]', ...
    [MarginDecomp.loc_shift]',[MarginDecomp.mean_loc_gap]', ...
    [MarginStats.d_median]',[MarginStats.d_iqr]', ...
    [ExtremeStats.d_cold]',[ExtremeStats.d_warm]',[ExtremeStats.d_total]', ...
    'VariableNames',{'Forcing','dMean','dLoc','dMeanMinusLoc','dMedian', ...
    'dIQR','dColdMass','dWarmMass','dExtreme'});

BenchmarkComparison = frstat_benchmark_comparison( ...
    y_grid,yvar0,xvar,est_all,Fnames);

coverageConstruction = repelem(["Baseline"; "Hemispheric equal weight"; ...
    "Fixed grid ("+string(fixedGridCells)+" cells)"],nForce);
coverageForcing = repmat(Fnames(:),3,1);
coverageMean = [[MarginStats.dmu_full]'; HemisphericCheck.table.dMean; ...
    FixedGridCheck.table.dMean];
coverageLocation = [[MarginDecomp.loc_shift]'; HemisphericCheck.table.dLoc; ...
    FixedGridCheck.table.dLoc];
coverageIQR = [[MarginStats.d_iqr]'; HemisphericCheck.table.dIQR; ...
    FixedGridCheck.table.dIQR];
coverageTail = [[ExtremeStats.d_warm]'; HemisphericCheck.table.dTailMass; ...
    FixedGridCheck.table.dTailMass];
coverageTotal = [[ExtremeStats.d_total]'; HemisphericCheck.table.dTotalExtreme; ...
    FixedGridCheck.table.dTotalExtreme];
coveragePValue = [repmat(BetweenCointMC.version2.pval,nForce,1); ...
    repmat(HemisphericCheck.between.version2.pval,nForce,1); ...
    repmat(FixedGridCheck.between.version2.pval,nForce,1)];
PaperCoverageTable = table(coverageConstruction,coverageForcing,coverageMean, ...
    coverageLocation,coverageIQR,coverageTail,coverageTotal,coveragePValue, ...
    'VariableNames',{'Construction','Forcing','dMean','dLoc','dIQR', ...
    'dWarmMass','dExtreme','BetweenPValue'});

%% Figure 1: Temperature-anomaly densities and anthropogenic forcing
Figure1Handle = frstat_show_figure1( ...
    y_grid,yvar0,expectedYears,xvar,Fnames);

%% Figure 2: Functional responses and implied densities
Figure2Handle = frstat_show_figure2( ...
    y_grid,yvar0,CIband,DensityCI,est_all,xvar,Fnames,ShockDeltaX);

%% Figure 3: Local probability-mass responses
Figure3Handle = frstat_show_figure3(LocalLossCurve,LocalLossCI,Fnames);

%% Table 1: Anthropogenic forcing portfolios
fprintf('\nTABLE 1. ANTHROPOGENIC FORCING PORTFOLIOS\n');
disp(ForcingAggregation);

%% Table 2: Margin responses
fprintf('\nTABLE 2. MARGIN RESPONSES\n');
disp(PaperMarginTable);

%% Table 3: Benchmark comparison
fprintf('\nTABLE 3. BENCHMARK COMPARISON\n');
disp(BenchmarkComparison);

%% Supplement Table 1: Coverage diagnostics
fprintf('\nSUPPLEMENT TABLE 1. COVERAGE DIAGNOSTICS\n');
disp(PaperCoverageTable);

%% Reported diagnostics
fprintf('\nREPORTED DIAGNOSTICS\n');
fprintf('F1/F2 level correlation: %.4f\n',F1F2Correlation);
fprintf('Within-forcing statistics: %.4f, %.4f; p-value: %.4f\n', ...
    WithinCointTest.coint_teststat(1),WithinCointTest.coint_teststat(2), ...
    WithinCointTest.pvalue);
fprintf('Between T_K: %.4f; 5%% critical value: %.4f; p-value: %.4f\n', ...
    BetweenCointMC.version1.teststat,BetweenCointMC.version1.qu, ...
    BetweenCointMC.version1.pval);
fprintf('Between T_V: %.4f; 5%% critical value: %.4f; p-value: %.4f\n', ...
    BetweenCointMC.version2.teststat,BetweenCointMC.version2.qu, ...
    BetweenCointMC.version2.pval);
fprintf('Common response: %.2f; 95%% critical value: %.2f; p-value: %.4f\n', ...
    CommonResponseTest.observedStatistic,CommonResponseTest.criticalValue95, ...
    CommonResponseTest.pValue);
fprintf('Retained CLR variation: %.1f%%; retained covariance trace: %.1f%%\n', ...
    100*CommonResponseSensitivity.responseSpaceShare, ...
    100*CommonResponseSensitivity.covarianceTraceShare);
fprintf('J-sensitivity p-value range: %.4f--%.4f\n', ...
    CommonResponseSensitivity.JPValueRange);
fprintf('M-sensitivity p-value range: %.4f--%.4f\n', ...
    CommonResponseSensitivity.MPValueRange);
fprintf('Observed cells: %d of %d in %d; %d in %d\n', ...
    observedCells(1),totalGridCells,coverageYears(1), ...
    observedCells(end),coverageYears(end));
fprintf('Fixed grid: %d cells (%.1f%%); median annual cell-months: %d and %d\n', ...
    fixedGridCells,100*fixedGridShare,medianObservedCellMonths, ...
    medianFixedGridCellMonths);
