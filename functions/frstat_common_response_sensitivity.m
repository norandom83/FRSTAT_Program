function out = frstat_common_response_sensitivity( ...
    yGrid,yvar0,xvar,S,nRep,seed,baselineTest)
%FRSTAT_COMMON_RESPONSE_SENSITIVITY Appendix common-response calculations.

jGrid = [8,12,16,20,24,28,32,40];
mGrid = [3,5,7,10,15,20];
fit20 = local_fit(yGrid,yvar0,xvar,20);

jPValue = NaN(size(jGrid));
jStatistic = NaN(size(jGrid));
for ii = 1:numel(jGrid)
    fit = local_fit(yGrid,yvar0,xvar,jGrid(ii));
    jStatistic(ii) = fit.observedStatistic;
    if jGrid(ii) == 20
        jPValue(ii) = baselineTest.pValue;
    else
        simulated = local_simulated_statistics(fit,7,S.bmgrid,nRep,seed);
        jPValue(ii) = (1+sum(simulated >= fit.observedStatistic))/(nRep+1);
    end
end

mPValue = NaN(size(mGrid));
for ii = 1:numel(mGrid)
    if mGrid(ii) == 7
        mPValue(ii) = baselineTest.pValue;
    else
        simulated = local_simulated_statistics( ...
            fit20,mGrid(ii),S.bmgrid,nRep,seed);
        mPValue(ii) = ...
            (1+sum(simulated >= fit20.observedStatistic))/(nRep+1);
    end
end

out = struct();
out.responseSpaceShare = fit20.retainedVariationShare;
out.covarianceTraceShare = sum(fit20.evals(1:7))/sum(fit20.evals);
out.bandwidth = fit20.bandwidth;
out.replications = nRep;
out.seed = seed;
out.J = table(jGrid(:),jStatistic(:),jPValue(:), ...
    'VariableNames',{'BasisDimension','Statistic','PValue'});
out.M = table(mGrid(:),repmat(fit20.observedStatistic,numel(mGrid),1), ...
    mPValue(:),'VariableNames',{'TruncationDimension','Statistic','PValue'});
out.JPValueRange = [min(jPValue),max(jPValue)];
out.MPValueRange = [min(mPValue),max(mPValue)];
end

function fit = local_fit(yGrid,yvar0,xvar,J)
nGrid = numel(yGrid);
weights = ones(nGrid,1)/(nGrid-1);
weights([1,end]) = weights([1,end])/2;
logDensity = log(yvar0);
clrY = logDensity-weights'*logDensity;
basis = local_fourier_basis(nGrid,J);
scores = basis'*(clrY.*weights);
scores = scores-mean(scores,2);

XX = xvar-mean(xvar,1);
betaPre = ((XX'*XX)\(XX'*scores'))';
residual = scores-betaPre*XX';
XD = diff(XX);
UD = residual(:,2:end)';
h = round(size(XD,1)^(1/4));

omegaXX2 = lr_var(XD,XD,2,h).omega;
omegaXX2 = (omegaXX2+omegaXX2')/2;
omegaXU2 = lr_var(UD,XD,2,h).omega;
omegaXU1 = lr_var(UD,XD,1,h).omega;
omegaXX1 = lr_var(XD,XD,1,h).omega;
upsilon = omegaXU1-omegaXU2/omegaXX2*omegaXX1;
XX2 = XX(2:end,:);
XX2 = XX2-mean(XX2,1);
nEst = size(XX2,1);
Z1 = scores(:,2:end)-omegaXU2/omegaXX2*XD';
beta = ((XX2'*XX2)\(XX2'*Z1'-nEst*upsilon'))';

omegaUU2 = lr_var(UD,UD,2,h).omega;
omegaUN = omegaUU2-omegaXU2/omegaXX2*omegaXU2';
omegaUN = (omegaUN+omegaUN')/2;
[vecUN,valUN] = eig(omegaUN,'vector');
[evals,indexUN] = sort(valUN,'descend');
vecUN = vecUN(:,indexUN);
[vecX,valX] = eig(omegaXX2,'vector');
[evalsx,indexX] = sort(valX,'descend');
vecX = vecX(:,indexX);

estimatedResponse = basis*beta;
estimatedDifference = estimatedResponse(:,1)-estimatedResponse(:,2);
fit = struct();
fit.yGrid = yGrid;
fit.evecs = basis*vecUN;
fit.evals = evals;
fit.evecsx = vecX;
fit.evalsx = evalsx;
fit.bandwidth = h;
fit.observedStatistic = ...
    nEst*sqrt(trapz(yGrid,estimatedDifference.^2));
fit.retainedVariationShare = ...
    sum(scores(2:end,:).^2,'all')*(yGrid(end)-yGrid(1))/ ...
    sum(trapz(yGrid,(clrY-mean(clrY,2)).^2,1));
end

function statistics = local_simulated_statistics(fit,M,bmgrid,nRep,seed)
priorRng = rng;
cleanupRng = onCleanup(@() rng(priorRng));
rng(seed,'twister');

responseVectors = fit.evecs(:,1:M);
responseSD = sqrt(max(fit.evals(1:M),0));
predictorMap = fit.evecsx*diag(sqrt(max(fit.evalsx,0)));
gram = NaN(M,M);
for ii = 1:M
    for jj = 1:M
        gram(ii,jj) = trapz( ...
            fit.yGrid,responseVectors(:,ii).*responseVectors(:,jj));
    end
end

statistics = NaN(nRep,1);
batchSize = 250;
dtRoot = sqrt(1/bmgrid);
for first = 1:batchSize:nRep
    last = min(first+batchSize-1,nRep);
    nBatch = last-first+1;
    innovations = randn(bmgrid,M+2,nBatch);
    responseIncrements = ...
        dtRoot*permute(innovations(:,1:M,:),[2,1,3]);
    responseIncrements = responseSD.*responseIncrements;

    predictorIncrements = ...
        dtRoot*permute(innovations(:,M+1:M+2,:),[2,1,3]);
    predictorMotion = ...
        cat(2,zeros(2,1,nBatch),cumsum(predictorIncrements,2));
    predictorMotion = predictorMotion(:,1:bmgrid,:);
    predictorLimit = pagemtimes(predictorMap,predictorMotion);
    predictorLimit = predictorLimit-mean(predictorLimit,2);

    denominator = ...
        pagemtimes(predictorLimit,'none',predictorLimit,'transpose')/bmgrid;
    aa = denominator(1,1,:);
    bb = denominator(1,2,:);
    cc = denominator(2,1,:);
    dd = denominator(2,2,:);
    determinant = aa.*dd-bb.*cc;
    v1 = (dd+bb)./determinant;
    v2 = -(cc+aa)./determinant;
    projectedPredictor = ...
        predictorLimit(1,:,:).*v1+predictorLimit(2,:,:).*v2;
    coefficientError = sum(responseIncrements.*projectedPredictor,2);
    gramError = pagemtimes(gram,coefficientError);
    statisticSquared = sum(coefficientError.*gramError,1);
    statistics(first:last) = ...
        reshape(sqrt(max(statisticSquared,0)),[],1);
end
end

function basis = local_fourier_basis(nGrid,J)
t = (0:(nGrid-1))'/(nGrid-1);
inner = @(f,g) trapz(t,f.*g);
raw = NaN(nGrid,J);
for ii = 1:(J/2)
    sine = sqrt(2)*sin(2*pi*ii*t);
    cosine = sqrt(2)*cos(2*pi*ii*t);
    raw(:,2*ii-1) = sine/sqrt(inner(sine,sine));
    raw(:,2*ii) = cosine/sqrt(inner(cosine,cosine));
end
basis = [ones(nGrid,1),raw];
basis(:,1) = basis(:,1)/sqrt(inner(basis(:,1),basis(:,1)));
for ii = 2:size(basis,2)
    for jj = 1:size(basis,2)
        if jj ~= ii
            basis(:,ii) = basis(:,ii)- ...
                inner(basis(:,ii),basis(:,jj))/ ...
                inner(basis(:,jj),basis(:,jj))*basis(:,jj);
        end
    end
end
for ii = 1:size(basis,2)
    basis(:,ii) = basis(:,ii)/sqrt(inner(basis(:,ii),basis(:,ii)));
end
end
