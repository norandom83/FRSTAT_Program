function out = frstat_reestimate_check(y_grid, yvar0, xvar, S, ...
    nHarmonics, seed, label)
%FRSTAT_REESTIMATE_CHECK Re-estimate the full specification on another density construction.

nGrid = numel(y_grid);
logDensity = log(yvar0);
trapWeights = local_trapezoid_weights(nGrid);
clrMean = trapWeights' * logDensity;
clrY = logDensity - clrMean;

basis = local_fourier_basis(nGrid,nHarmonics);
scores = basis' * (clrY .* trapWeights);
scores = scores - mean(scores,2);
XX = xvar - mean(xvar,1);

betaPre = ((XX'*XX) \ (XX'*scores'))';
residual = scores-betaPre*XX';
XD = diff(XX);
UD = residual(:,2:end)';
kerband = round(size(XD,1)^(1/4));

omegaXX2 = lr_var(XD,XD,2,kerband).omega;
omegaXX2 = (omegaXX2+omegaXX2')/2;
if rcond(omegaXX2)<1e-12
    error('The forcing long-run covariance is numerically singular.');
end
omegaXU2 = lr_var(UD,XD,2,kerband).omega;
omegaXU1 = lr_var(UD,XD,1,kerband).omega;
omegaXX1 = lr_var(XD,XD,1,kerband).omega;
upsilon = omegaXU1-omegaXU2/omegaXX2*omegaXX1;
Z1 = scores(:,2:end)-omegaXU2/omegaXX2*XD';
XX2 = XX(2:end,:);
XX2 = XX2-mean(XX2,1);
nEst = size(XX2,1);
beta = ((XX2'*XX2) \ (XX2'*Z1'-nEst*upsilon'))';

nForce = size(xvar,2);
deltaX = std(xvar,0,1)';
fRef = build_reference_density_seed(y_grid,yvar0);
dMean=NaN(nForce,1); dLoc=NaN(nForce,1); dMedian=NaN(nForce,1);
dIQR=NaN(nForce,1); dTail=NaN(nForce,1); dTotal=NaN(nForce,1);
for ip=1:nForce
    est=(beta(:,ip)'*basis')';
    [ms,md]=frstat_margin_step(y_grid,yvar0,est,deltaX(ip),S,"F"+string(ip));
    ex=frstat_extreme_step(y_grid,fRef,est,deltaX(ip),md.loc_shift,S,"F"+string(ip));
    dMean(ip)=ms.dmu_full; dLoc(ip)=md.loc_shift; dMedian(ip)=ms.d_median;
    dIQR(ip)=ms.d_iqr; dTail(ip)=ms.d_tail; dTotal(ip)=ex.d_total;
end

between=frstat_between_coint_mc_tests(residual,xvar,'Alpha',0.05, ...
    'Bmc',S.testReplications,'McGrid',S.testMcGrid,'TestDemean',true,'Seed',seed);
out=struct('label',string(label),'support',[y_grid(1),y_grid(end)], ...
    'y_grid',y_grid,'yvar0',yvar0,'beta',beta,'between',between);
out.table=table("F"+string((1:nForce)'),dMean,dLoc,dMedian,dIQR,dTail,dTotal, ...
    'VariableNames',{'Forcing','dMean','dLoc','dMedian','dIQR','dTailMass','dTotalExtreme'});
end

function weights = local_trapezoid_weights(nGrid)
weights = ones(nGrid,1)/(nGrid-1);
weights([1,end]) = weights([1,end])/2;
end

function basis = local_fourier_basis(nGrid,nHarmonics)
t=(0:(nGrid-1))'/(nGrid-1);
inner=@(f,g) trapz(t,f.*g);
raw=NaN(nGrid,nHarmonics);
for i=1:(nHarmonics/2)
    s=sqrt(2)*sin(2*pi*i*t); c=sqrt(2)*cos(2*pi*i*t);
    raw(:,2*i-1)=s/sqrt(inner(s,s)); raw(:,2*i)=c/sqrt(inner(c,c));
end
basis=[ones(nGrid,1),raw];
basis(:,1)=basis(:,1)/sqrt(inner(basis(:,1),basis(:,1)));
for i=2:size(basis,2)
    for j=1:size(basis,2)
        if j~=i
            basis(:,i)=basis(:,i)-(inner(basis(:,i),basis(:,j))/inner(basis(:,j),basis(:,j)))*basis(:,j);
        end
    end
end
for i=1:size(basis,2)
    basis(:,i)=basis(:,i)/sqrt(inner(basis(:,i),basis(:,i)));
end
end
