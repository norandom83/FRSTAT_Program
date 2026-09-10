function FigureHandle = frstat_show_figure2(y_grid, yvar0, CIband, ...
    DensityCI, est_all, xvar, displayNames, shockInput)
%% Display Figure 2: CLR and implied-density responses.

y_grid = y_grid(:);
nForce = size(est_all, 2);
if size(est_all, 1) ~= numel(y_grid) || numel(CIband) < nForce || ...
        numel(DensityCI) < nForce
    error('Figure 2 response inputs do not match the temperature grid.');
end
if size(xvar, 2) < nForce
    error('Figure 2 requires one forcing series for each response function.');
end

if nargin < 7 || isempty(displayNames)
    displayNames = "F" + string(1:nForce);
else
    displayNames = string(displayNames(:))';
end
if numel(displayNames) ~= nForce
    error('The number of forcing names must match the response functions.');
end

if nargin < 8 || isempty(shockInput)
    shockDeltaX = std(xvar(:,1:nForce), 0, 1)';
elseif isstruct(shockInput) && isfield(shockInput, 'deltaX')
    shockDeltaX = shockInput.deltaX(:);
else
    shockDeltaX = shockInput(:);
end
if numel(shockDeltaX) ~= nForce || any(~isfinite(shockDeltaX))
    error('Figure 2 requires one finite shock magnitude for each forcing.');
end

forcingColor = forcing_colors(nForce);
p_ref = build_reference_density_seed(y_grid, yvar0);

responseEstimate = cell(nForce, 1);
responseLower = cell(nForce, 1);
responseUpper = cell(nForce, 1);
counterfactualDensity = NaN(numel(y_grid), nForce);
densityLower = NaN(numel(y_grid), nForce);
densityUpper = NaN(numel(y_grid), nForce);
allResponseValues = [];

for ip = 1:nForce
    requiredField = {'estimate', 'lquan', 'uquan'};
    if ~all(isfield(CIband(ip), requiredField))
        error('Each CIband element must contain estimate, lquan, and uquan.');
    end
    estimate = CIband(ip).estimate(:);
    lowerError = CIband(ip).lquan(:);
    upperError = CIband(ip).uquan(:);
    if any([numel(estimate), numel(lowerError), numel(upperError)] ~= numel(y_grid))
        error('Each Figure 2 response band must match the temperature grid.');
    end

    deltaX = shockDeltaX(ip);
    responseEstimate{ip} = estimate * deltaX;
    responseLower{ip} = min((estimate + lowerError) * deltaX, ...
        (estimate + upperError) * deltaX);
    responseUpper{ip} = max((estimate + lowerError) * deltaX, ...
        (estimate + upperError) * deltaX);
    allResponseValues = [allResponseValues; responseEstimate{ip}; ...
        responseLower{ip}; responseUpper{ip}]; %#ok<AGROW>

    densityFields = {'estimate', 'lower', 'upper'};
    if ~all(isfield(DensityCI(ip), densityFields))
        error('Each DensityCI element must contain estimate, lower, and upper.');
    end
    counterfactualDensity(:,ip) = DensityCI(ip).estimate(:);
    densityLower(:,ip) = DensityCI(ip).lower(:);
    densityUpper(:,ip) = DensityCI(ip).upper(:);
    if any([numel(counterfactualDensity(:,ip)), ...
            numel(densityLower(:,ip)), numel(densityUpper(:,ip))] ~= ...
            numel(y_grid)) || any(~isfinite([counterfactualDensity(:,ip); ...
            densityLower(:,ip); densityUpper(:,ip)])) || ...
            any(densityLower(:,ip) > densityUpper(:,ip))
        error('Each density interval must be finite and match the temperature grid.');
    end
end

commonResponseLimit = padded_limits(allResponseValues, 0.08, true);
commonDensityMaximum = ceil(1.06*max([p_ref; densityUpper(:)])/0.05)*0.05;
commonDensityMaximum = max(commonDensityMaximum, 0.55);

FigureHandle = struct();
FigureHandle.CLR = gobjects(nForce, 1);
FigureHandle.Density = gobjects(nForce, 1);
FigureHandle.CommonDensityYLim = [0 commonDensityMaximum];

for ip = 1:nForce
    if shockDeltaX(ip) >= 0
        shockLabel = '+1 s.d. change';
    else
        shockLabel = '-1 s.d. change';
    end

    %% CLR response panel
    FigureHandle.CLR(ip) = figure('Color', 'w', ...
        'Position', [100 100 1000 720]);
    hold on;
    hEstimate = plot(y_grid, responseEstimate{ip}, 'k-', 'LineWidth', 3.0);
    hLimit = plot(y_grid, responseLower{ip}, 'k--', 'LineWidth', 2.5);
    plot(y_grid, responseUpper{ip}, 'k--', 'LineWidth', 2.5, ...
        'HandleVisibility', 'off');
    xlabel('Temperature anomaly (\circC)', 'fontsize', 27, 'fontweight', 'b');
    ylabel('CLR perturbation', 'fontsize', 27, 'fontweight', 'b');
    title("CLR response: " + displayNames(ip) + " (" + shockLabel + ")", ...
        'fontsize', 24, 'fontweight', 'b');
    legend([hEstimate hLimit], ...
        {'Local-average projection', 'Marginal 95% limits'}, ...
        'Location', 'best', 'box', 'off', 'FontSize', 18);
    hold off; grid on;
    xlim([y_grid(1) y_grid(end)]);
    if ip == 2
        ylim([-1 0.5]);
    else
        ylim(commonResponseLimit);
    end
    set(gca, 'LineWidth', 2.5, 'box', 'on', 'TickLength', [0 0], ...
        'FontSize', 24, 'fontweight', 'b');

    %% Implied-density response panel
    FigureHandle.Density(ip) = figure('Color', 'w', ...
        'Position', [100 100 1150 535]);
    ax = axes(FigureHandle.Density(ip), ...
        'Position', [0.087 0.15 0.889 0.765]);
    hold(ax, 'on');
    fill(ax, [y_grid; flipud(y_grid)], ...
        [densityLower(:,ip); flipud(densityUpper(:,ip))], ...
        forcingColor(ip,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15, ...
        'HandleVisibility', 'off');
    hDensityLimit = plot(ax, y_grid, densityLower(:,ip), ':', ...
        'Color', forcingColor(ip,:), 'LineWidth', 1.7);
    plot(ax, y_grid, densityUpper(:,ip), ':', ...
        'Color', forcingColor(ip,:), 'LineWidth', 1.7, ...
        'HandleVisibility', 'off');
    hReference = plot(ax, y_grid, p_ref, '-', ...
        'Color', [0.27 0.27 0.27], 'LineWidth', 2.3);
    hCounterfactual = plot(ax, y_grid, counterfactualDensity(:,ip), '-', ...
        'Color', forcingColor(ip,:), 'LineWidth', 2.8);
    xlabel(ax, 'Temperature anomaly (\circC)', ...
        'FontSize', 19, 'FontWeight', 'bold');
    ylabel(ax, 'Density', 'FontSize', 19, 'FontWeight', 'bold');
    title(ax, "Density response: " + displayNames(ip), ...
        'FontSize', 21, 'FontWeight', 'bold');
    legend(ax, [hReference hCounterfactual hDensityLimit], ...
        {'Reference density', char(shockLabel + " to " + displayNames(ip)), ...
         'Pointwise 95% CI'}, ...
        'Location', 'northwest', 'Box', 'off', 'FontSize', 15, ...
        'AutoUpdate', 'off');
    xlim(ax, [y_grid(1) y_grid(end)]);
    ylim(ax, [0 commonDensityMaximum]);
    xticks(ax, -5:1:5);
    yticks(ax, 0:0.1:commonDensityMaximum);
    set(ax, 'FontName', 'Arial', 'FontSize', 17, 'LineWidth', 1.1, ...
        'Box', 'on', 'TickDir', 'out', 'TickLength', [0.004 0.004], ...
        'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12, 'Layer', 'top');
    hold(ax, 'off');
end
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
