function FigureHandle = frstat_show_figure3(LocalLossCurve, LocalLossCI, displayNames)
%% Display Figure 3: local probability-mass decompositions.

nForce = numel(LocalLossCurve);
if nargin < 3 || isempty(displayNames)
    displayNames = strings(1, nForce);
    for ip = 1:nForce
        if isfield(LocalLossCurve(ip), 'forcing') && ...
                strlength(string(LocalLossCurve(ip).forcing)) > 0
            displayNames(ip) = string(LocalLossCurve(ip).forcing);
        else
            displayNames(ip) = "F" + string(ip);
        end
    end
else
    displayNames = string(displayNames(:))';
end
if numel(displayNames) ~= nForce
    error('The number of forcing names must match the local-response profiles.');
end

forcingColor = forcing_colors(nForce);
xMinimum = inf;
xMaximum = -inf;
allLocalValues = [];

for ip = 1:nForce
    requiredField = {'r_grid', 'dLoss_full', 'dLoss_loc', 'dLoss_dist'};
    if ~all(isfield(LocalLossCurve(ip), requiredField))
        error('Each local-response profile must contain its grid and three components.');
    end
    rGrid = LocalLossCurve(ip).r_grid(:);
    fullResponse = LocalLossCurve(ip).dLoss_full(:);
    locationResponse = LocalLossCurve(ip).dLoss_loc(:);
    distributionResponse = LocalLossCurve(ip).dLoss_dist(:);
    if any([numel(fullResponse), numel(locationResponse), ...
            numel(distributionResponse)] ~= numel(rGrid))
        error('Each Figure 3 response component must match its local grid.');
    end

    xMinimum = min(xMinimum, min(rGrid));
    xMaximum = max(xMaximum, max(rGrid));
    allLocalValues = [allLocalValues; fullResponse; ...
        locationResponse; distributionResponse]; %#ok<AGROW>
    [hasBand, ~, bandLower, bandUpper] = local_loss_band(LocalLossCI, ip);
    if hasBand
        allLocalValues = [allLocalValues; bandLower; bandUpper]; %#ok<AGROW>
    end
end

commonYLim = padded_limits(allLocalValues, 0.12, true);
FigureHandle = struct();
FigureHandle.Panel = gobjects(nForce, 1);
FigureHandle.CommonXLim = [xMinimum xMaximum];
FigureHandle.CommonYLim = commonYLim;

for ip = 1:nForce
    rGrid = LocalLossCurve(ip).r_grid(:);
    FigureHandle.Panel(ip) = figure('Color', 'w', ...
        'Position', [100 100 1400 650]);
    hold on;

    [hasBand, bandGrid, bandLower, bandUpper] = ...
        local_loss_band(LocalLossCI, ip);
    if hasBand
        hBand = patch([bandGrid; flipud(bandGrid)], ...
            [bandLower; flipud(bandUpper)], [0.65 0.65 0.65], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.35);
    else
        hBand = patch(NaN, NaN, [0.65 0.65 0.65], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.35);
    end

    hFull = plot(rGrid, LocalLossCurve(ip).dLoss_full(:), '-', ...
        'LineWidth', 4.0, 'Color', forcingColor(ip,:));
    hLocation = plot(rGrid, LocalLossCurve(ip).dLoss_loc(:), '--', ...
        'LineWidth', 3.5, 'Color', [0.15 0.15 0.15]);
    hDistribution = plot(rGrid, LocalLossCurve(ip).dLoss_dist(:), '-.', ...
        'LineWidth', 3.5, 'Color', [0.45 0.45 0.45]);
    yline(0, 'k-', 'LineWidth', 1.0, 'HandleVisibility', 'off');

    xlabel('Temperature anomaly (\circC)', 'fontsize', 28, 'fontweight', 'b');
    ylabel('\Delta probability mass', 'fontsize', 27, 'fontweight', 'b');
    title("Local probability-mass change: " + displayNames(ip), ...
        'fontsize', 24, 'fontweight', 'b');
    legend([hFull hLocation hDistribution hBand], ...
        {'\Delta(full)', '\Delta(loc-only)', '\Delta(dist)', ...
         '95% pointwise interval: full'}, ...
        'Location', 'best', 'box', 'off', 'FontSize', 19);
    hold off; grid on;
    xlim([xMinimum xMaximum]);
    ylim(commonYLim);
    ax = gca;
    ax.YAxis.Exponent = 0;
    set(gca, 'LineWidth', 2.5, 'box', 'on', 'TickLength', [0 0], ...
        'FontSize', 24, 'fontweight', 'b');
end
end

function [hasBand, rGrid, lowerBand, upperBand] = local_loss_band(LocalLossCI, ip)
hasBand = false;
rGrid = [];
lowerBand = [];
upperBand = [];

if ~isstruct(LocalLossCI) || isempty(LocalLossCI) || ...
        ~isfield(LocalLossCI, 'r_grid') || ~isfield(LocalLossCI, 'bands') || ...
        numel(LocalLossCI.bands) < ip
    return;
end

band = LocalLossCI.bands(ip);
if ~isstruct(band) || ~isfield(band, 'lower') || ~isfield(band, 'upper')
    return;
end

rGrid = LocalLossCI.r_grid(:);
lowerBand = band.lower(:);
upperBand = band.upper(:);
nObservation = min([numel(rGrid), numel(lowerBand), numel(upperBand)]);
if nObservation < 2
    rGrid = [];
    lowerBand = [];
    upperBand = [];
    return;
end

rGrid = rGrid(1:nObservation);
lowerBand = lowerBand(1:nObservation);
upperBand = upperBand(1:nObservation);
valid = isfinite(rGrid) & isfinite(lowerBand) & isfinite(upperBand);
if nnz(valid) < 2
    rGrid = [];
    lowerBand = [];
    upperBand = [];
    return;
end

rGrid = rGrid(valid);
lowerValue = min(lowerBand(valid), upperBand(valid));
upperValue = max(lowerBand(valid), upperBand(valid));
[rGrid, sortIndex] = sort(rGrid);
lowerBand = lowerValue(sortIndex);
upperBand = upperValue(sortIndex);
hasBand = true;
end

function plotColor = forcing_colors(nForce)
baseColor = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; ...
    0.4660 0.6740 0.1880; 0.4940 0.1840 0.5560; ...
    0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];
plotColor = lines(nForce);
plotColor(1:min(nForce, size(baseColor,1)), :) = ...
    baseColor(1:min(nForce, size(baseColor,1)), :);
end

function lim = padded_limits(values, padFraction, includeZero)
values = values(isfinite(values));
if isempty(values)
    lim = [-1 1];
    return;
end

lowerValue = min(values);
upperValue = max(values);
if includeZero
    lowerValue = min(lowerValue, 0);
    upperValue = max(upperValue, 0);
end

valueRange = upperValue - lowerValue;
if valueRange <= 0
    valueRange = max(abs([lowerValue upperValue 1]));
end
padding = padFraction * valueRange;
lim = [lowerValue - padding, upperValue + padding];
end
