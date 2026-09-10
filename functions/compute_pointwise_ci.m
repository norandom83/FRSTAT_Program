function [lquan, uquan, projectionMatrix] = compute_pointwise_ci(bounds, y_grid, stepsize)
    G = length(y_grid);
    III = zeros(G, G);

    for i = 1:G
        II = zeros(1, G);
        lo = max(i-stepsize,1);
        hi = min(i+stepsize,G);
        windowWeights = ones(hi-lo+1,1);
        if numel(windowWeights) > 1
            windowWeights([1,end]) = 0.5;
        end
        II(lo:hi) = (windowWeights / sum(windowWeights))';
        III(i,:) = II / sum(II);
    end

    ciq   = zeros(size(III,1), size(bounds,2));
    lquan = zeros(size(III,1),1);
    uquan = zeros(size(III,1),1);

    for j = 1:size(III,1)
        for i = 1:size(bounds,2)
            ciq(j,i) = sum(III(j,:)' .* bounds(:,i));
        end
        qLow = quantile(ciq(j,:), 0.025);
        qHigh = quantile(ciq(j,:), 0.975);
        % bounds simulate estimation error, so invert its quantiles.
        lquan(j) = -qHigh;
        uquan(j) = -qLow;
    end
    projectionMatrix = III;
end
