function FigureHandle = frstat_show_figure1(y_grid, yvar0, years, ...
    xvar, displayNames)
%% Display Figure 1: temperature densities and anthropogenic forcing.

y_grid = y_grid(:);
T = size(yvar0, 2);

if nargin < 3 || isempty(years)
    years = (1:T)';
else
    years = years(:);
end

if size(yvar0, 1) ~= numel(y_grid) || numel(years) ~= T || ...
        size(xvar, 1) ~= T
    error('Figure 1 inputs must have a common temperature grid and time dimension.');
end
nForce = size(xvar, 2);
if nargin < 5 || isempty(displayNames)
    displayNames = "F" + string(1:nForce);
else
    displayNames = string(displayNames(:))';
end
if numel(displayNames) ~= nForce
    error('The number of forcing names must match the number of forcing series.');
end

FigureHandle = struct();

%% Distribution surface
FigureHandle.DistributionSurface = figure('Color', 'w', ...
    'Position', [100 100 1200 800]);
[Xsurf, Ysurf] = meshgrid(y_grid, years);
surf(Xsurf, Ysurf, yvar0', 'EdgeColor', 'none');
view(135, 30);
colormap(parula);
colorbar;
xlabel('Temperature anomaly (\circC)', 'fontsize', 28, 'fontweight', 'b');
ylabel('Year', 'fontsize', 28, 'fontweight', 'b');
zlabel('Density', 'fontsize', 28, 'fontweight', 'b');
title('Evolution of the temperature anomaly distribution', ...
    'fontsize', 28, 'fontweight', 'b');
axis tight;
set(gca, 'LineWidth', 2.5, 'box', 'on', 'TickLength', [0 0], ...
    'FontSize', 25, 'fontweight', 'b');

%% Anthropogenic forcing portfolios
forcingColor = forcing_colors(nForce);
FigureHandle.Forcing = figure('Color', 'w', 'Position', [100 100 1350 700]);
hold on;
lineHandle = gobjects(nForce, 1);
for jj = 1:nForce
    lineHandle(jj) = plot(years, xvar(:,jj), '-', ...
        'Color', forcingColor(jj,:), 'LineWidth', 4.0);
end
legend(lineHandle, cellstr(displayNames), 'Location', 'northwest', ...
    'box', 'off', 'FontSize', 22);
xlabel('Year', 'fontsize', 28, 'fontweight', 'b');
ylabel('Effective radiative forcing (W/m^2)', ...
    'fontsize', 28, 'fontweight', 'b');
title('Anthropogenic forcing portfolios', ...
    'fontsize', 28, 'fontweight', 'b');
hold off; grid on; axis tight;
set(gca, 'LineWidth', 2.5, 'box', 'on', 'TickLength', [0 0], ...
    'FontSize', 25, 'fontweight', 'b');

%% Selected annual densities
idx = unique([1, round(T/2), T], 'stable');
plotColor = [0 0.4470 0.7410; 0.20 0.20 0.20; 0.8500 0.3250 0.0980];
plotColor = plotColor(1:numel(idx), :);
legendText = cellstr(string(years(idx)));

FigureHandle.SelectedDensities = figure('Color', 'w', ...
    'Position', [100 100 1050 720]);
hold on;
for jj = 1:numel(idx)
    plot(y_grid, yvar0(:,idx(jj)), '-', 'Color', plotColor(jj,:), ...
        'LineWidth', 4.0);
end
xlabel('Temperature anomaly (\circC)', 'fontsize', 28, 'fontweight', 'b');
ylabel('Density', 'fontsize', 28, 'fontweight', 'b');
title('Global temperature anomaly distributions', ...
    'fontsize', 28, 'fontweight', 'b');
legend(legendText, 'Location', 'northwest', 'box', 'off', 'FontSize', 22);
hold off; grid on; axis tight;
set(gca, 'LineWidth', 2.5, 'box', 'on', 'TickLength', [0 0], ...
    'FontSize', 25, 'fontweight', 'b');

%% Mean-aligned annual densities
FigureHandle.MeanAlignedDensities = figure('Color', 'w', ...
    'Position', [100 100 1050 720]);
hold on;
for jj = 1:numel(idx)
    density = normalize_density_trapz(max(yvar0(:,idx(jj)), eps), y_grid);
    meanValue = trapz(y_grid, y_grid .* density);
    plot(y_grid - meanValue, density, '-', 'Color', plotColor(jj,:), ...
        'LineWidth', 4.0);
end
xlabel('Mean-aligned temperature anomaly (\circC)', ...
    'fontsize', 28, 'fontweight', 'b');
ylabel('Density', 'fontsize', 28, 'fontweight', 'b');
title('Mean-aligned anomaly distributions', ...
    'fontsize', 28, 'fontweight', 'b');
legend(legendText, 'Location', 'northwest', 'box', 'off', 'FontSize', 22);
hold off; grid on; axis tight;
set(gca, 'LineWidth', 2.5, 'box', 'on', 'TickLength', [0 0], ...
    'FontSize', 25, 'fontweight', 'b');

%% Warm-side probability reallocation
nBlock = round(T/3);
firstDensity = mean(yvar0(:,1:nBlock), 2, 'omitnan');
lastDensity = mean(yvar0(:,(2*nBlock+1):T), 2, 'omitnan');
firstDensity = normalize_density_trapz(max(firstDensity, eps), y_grid);
lastDensity = normalize_density_trapz(max(lastDensity, eps), y_grid);
densityChange = lastDensity - firstDensity;
firstCDF = cumtrapz(y_grid, firstDensity);
firstCDF = firstCDF ./ firstCDF(end);
upperTailStart = interp1(firstCDF, y_grid, 0.90, 'linear', 'extrap');

FigureHandle.WarmSideReallocation = figure('Color', 'w', ...
    'Position', [100 100 1050 720]);
hold on;
tailIndex = y_grid >= upperTailStart;
if any(tailIndex)
    yLimit = padded_limits(densityChange, 0.08, true);
    xTail = y_grid(tailIndex);
    hTail = patch([xTail(1); xTail(end); xTail(end); xTail(1)], ...
        [yLimit(1); yLimit(1); yLimit(2); yLimit(2)], [1.00 0.80 0.80], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.35);
else
    hTail = plot(NaN, NaN, '-', 'Color', [1.00 0.80 0.80], 'LineWidth', 8);
end
hChange = plot(y_grid, densityChange, 'k-', 'LineWidth', 4.0);
hZero = yline(0, 'k--', 'LineWidth', 1.8);
xlabel('Temperature anomaly (\circC)', 'fontsize', 28, 'fontweight', 'b');
ylabel('Density change', 'fontsize', 28, 'fontweight', 'b');
title('Warm-side probability reallocation', ...
    'fontsize', 28, 'fontweight', 'b');
legend([hChange hZero hTail], ...
    {'Last third - first third', 'Zero line', 'Upper-tail region'}, ...
    'Location', 'northwest', 'box', 'off', 'FontSize', 20);
hold off; grid on;
xlim([y_grid(1) y_grid(end)]);
ylim(padded_limits(densityChange, 0.08, true));
set(gca, 'LineWidth', 2.5, 'box', 'on', 'TickLength', [0 0], ...
    'FontSize', 25, 'fontweight', 'b');
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
