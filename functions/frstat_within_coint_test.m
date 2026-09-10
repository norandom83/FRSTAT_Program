function out = frstat_within_coint_test(xvar, CVx, alpha)
%FRSTAT_WITHIN_COINT_TEST Breitung diagnostic for the forcing vector.

XX = xvar - mean(xvar,1);
XX = XX(2:end,:);
T = size(XX,1);
nForce = size(XX,2);

BC1 = XX'*XX;
partialSums = cumsum(XX,1);
BC2 = partialSums'*partialSums;
eigenvalues = sort(real(eig(BC1/BC2)),'ascend');
statistics = cumsum(eigenvalues*T^2);
testStatistic = statistics(nForce);

criticalIndex = find(testStatistic <= CVx(nForce,:),1,'first');
if isempty(criticalIndex)
    pValue = 0;
else
    pValue = 1-criticalIndex/size(CVx,2);
end

out = struct('coint_teststat',statistics,'teststat_used',testStatistic, ...
    'pvalue',pValue,'present_at_alpha',pValue < alpha);
end
