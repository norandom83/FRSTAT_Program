function BenchmarkComparison = frstat_benchmark_comparison( ...
    y_grid, yvar0, xvar, est_all, Fnames)
%FRSTAT_BENCHMARK_COMPARISON Construct the scalar benchmark comparison.

y_grid = y_grid(:);
T = size(xvar,1);
nForce = size(xvar,2);
aggregateForcing = sum(xvar,2);

densityMean = NaN(T,1);
for t = 1:T
    density = max(yvar0(:,t),eps);
    density = density/trapz(y_grid,density);
    densityMean(t) = trapz(y_grid,y_grid.*density);
end

scalarFM = frstat_fm_mean_slopes(densityMean,aggregateForcing);
dXAggregate = std(aggregateForcing,0,1);
dMeanScalar = scalarFM.slope*dXAggregate;

% The joint vector change sums to one standard deviation of F1+F2.
dXVector = std(xvar,0,1)';
dXVector = dXAggregate*dXVector/sum(dXVector);

vectorFM = frstat_fm_mean_slopes(densityMean,xvar);
dMeanVector = vectorFM.slope(:).*dXVector;
dMeanVectorJoint = sum(dMeanVector);

referenceDensity = max(mean(yvar0,2),eps);
referenceDensity = referenceDensity/trapz(y_grid,referenceDensity);
referenceCLR = log(referenceDensity);
referenceCLR = referenceCLR - ...
    trapz(y_grid,referenceCLR)/(y_grid(end)-y_grid(1));
[mean0,IQR0,tail0,r95] = local_density_margins(y_grid,referenceDensity);

dMeanDensity = NaN(nForce,1);
dIQRDensity = NaN(nForce,1);
dTailDensity = NaN(nForce,1);
for ip = 1:nForce
    clrNew = referenceCLR + dXVector(ip)*est_all(:,ip);
    densityNew = exp(clrNew-max(clrNew));
    densityNew = densityNew/trapz(y_grid,densityNew);
    [meanNew,IQRNew,tailNew] = local_density_margins( ...
        y_grid,densityNew,r95);
    dMeanDensity(ip) = meanNew-mean0;
    dIQRDensity(ip) = IQRNew-IQR0;
    dTailDensity(ip) = tailNew-tail0;
end

clrJoint = referenceCLR + est_all*dXVector;
densityJoint = exp(clrJoint-max(clrJoint));
densityJoint = densityJoint/trapz(y_grid,densityJoint);
[meanJoint,IQRJoint,tailJoint] = local_density_margins( ...
    y_grid,densityJoint,r95);

Model = ["Aggregate scalar"; repmat("Vector-to-mean",nForce,1); ...
    "Vector-to-mean"; repmat("Vector-to-density",nForce,1); ...
    "Vector-to-density"];
Channel = ["F1+F2"; "Matched "+Fnames(:)+" component"; ...
    "Matched joint change"; "Matched "+Fnames(:)+" component"; ...
    "Matched joint change"];
DeltaMean = [dMeanScalar; dMeanVector; dMeanVectorJoint; ...
    dMeanDensity; meanJoint-mean0];
DeltaIQR = [NaN; NaN(nForce,1); NaN; dIQRDensity; IQRJoint-IQR0];
DeltaWarmMass = [NaN; NaN(nForce,1); NaN; dTailDensity; tailJoint-tail0];

BenchmarkComparison = table(Model,Channel,DeltaMean,DeltaIQR,DeltaWarmMass);
end

function [densityMean,IQR,tailMass,r95] = local_density_margins( ...
    y_grid,density,r95)

density = max(density(:),eps);
density = density/trapz(y_grid,density);
densityMean = trapz(y_grid,y_grid.*density);

cdf = cumtrapz(y_grid,density);
cdf = cdf/cdf(end);
[cdfUnique,index] = unique(cdf,'stable');
q25 = interp1(cdfUnique,y_grid(index),0.25,'linear','extrap');
q75 = interp1(cdfUnique,y_grid(index),0.75,'linear','extrap');
IQR = q75-q25;

if nargin < 3
    r95 = interp1(cdfUnique,y_grid(index),0.95,'linear','extrap');
end
tailMass = density_interval_mass(y_grid,density,r95,Inf);
end
