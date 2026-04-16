function [RL, RL1, RL2, RL3, RL4, RL5] = Return_level_plot_3(params, AMS, YEAR_Range, best_fit, Label, no, no1, no2, no3, no4)
%==========================================================================
% Return_level_plot_3
% Plots 6 return periods (2, 5, 10, 25, 50, 100 years) in non-stationary mode
% with solid AMS data line color-coded by tile index (first 2, next 2, last 2).
%
% Inputs:
%   params     : cell array of sampled parameters
%   AMS        : annual maxima series
%   YEAR_Range : vector of years
%   best_fit   : 'gev', 'gp', 'gamma', or 'lognormal'
%   Label      : y-axis label
%   no         : show CI bands (1=yes, 0=no)
%   no1, no2, no3, no4 : control title, legend, xlabel, ylabel
%
% Outputs:
%   RL, RL1, RL2, RL3, RL4, RL5 : Return level series for each return period
%==========================================================================


RP_COLORS = [...
    0.95 0.60 0.00;  % Orange - 2-year
    0.55 0.00 0.90;  % Purple - 5-year
    0.00 0.75 0.75;  % Teal   - 10-year
    0.80 0.20 0.60;  % Magenta - 25-year
    0.60 0.40 0.00;  % Brown  - 50-year
    0.50 0.00 0.13]; % Black  - 100-year

Return_Period = [2; 5; 10; 25; 50; 100];
X_fit = YEAR_Range;
RL_all = cell(1, length(Return_Period));



% ---- AMS Data Line Color by Tile ----
if ismember(no1, [1, 2])
    dataColor = [0.2, 0.5, 0.9];    % blue (tiles 1–2)
elseif ismember(no1, [3, 4])
    dataColor = [0.2, 0.5, 0.9];    % Blue (tiles 3–4)
elseif ismember(no1, [5, 6])
    dataColor = [0.9, 0.3, 0.3];    % Red (tiles 5–6)
else
    % --- Handle unexpected values gracefully ---
    disp(['⚠️ Warning: Unexpected no1 value = ', num2str(no1)]);
    dataColor = [0.2, 0.5, 0.9];    % Default to Blue for safety
end



% ---- Figure setup ----
set(gcf, 'Color', 'white');
set(gca, 'Color', 'white');
nexttile
hold on

% ---- Distribution Handling ----
switch lower(best_fit)
    case 'gev'
        shape = params{1}; scale_t = params{2}; mu_t = params{3};
        X1 = mean(shape); X2 = mean(scale_t); X3 = mean(mu_t);
        X1_95 = prctile(shape,95); X1_05 = prctile(shape,5);
        X2_95 = prctile(scale_t,95); X2_05 = prctile(scale_t,5);
        X3_95 = prctile(mu_t,95); X3_05 = prctile(mu_t,5);

        for i = 1:length(Return_Period)
            RLtemp = gevinv(1 - 1/Return_Period(i), X1, exp(X2), X3);
            if no == 1
                CI_05 = gevinv(1 - 1/Return_Period(i), X1_05, exp(X2_05), X3_05);
                CI_95 = gevinv(1 - 1/Return_Period(i), X1_95, exp(X2_95), X3_95);
                % optional CI shading
            end
            p(i) = plot(X_fit, RLtemp, 'LineStyle','-.','Color', RP_COLORS(i,:), 'LineWidth', 1.4);
            RL_all{i} = RLtemp;
        end

    case 'gp'
        shape = params{1}; scale_t = params{2};
        X1 = mean(shape); X2 = mean(scale_t);
        X1_95 = prctile(shape,95); X1_05 = prctile(shape,5);
        X2_95 = prctile(scale_t,95); X2_05 = prctile(scale_t,5);

        for i = 1:length(Return_Period)
            RLtemp = gpinv(1 - 1/Return_Period(i), X1, exp(X2));
            if no == 1
                CI_05 = gpinv(1 - 1/Return_Period(i), X1_05, exp(X2_05));
                CI_95 = gpinv(1 - 1/Return_Period(i), X1_95, exp(X2_95));
            end
            p(i) = plot(X_fit, RLtemp, '-.', 'Color', RP_COLORS(i,:), 'LineWidth', 1.4);
            RL_all{i} = RLtemp;
        end

    case 'gamma'
        shape = params{1}; scale_t = params{2};
        X1 = mean(shape); X2 = mean(scale_t);
        X1_95 = prctile(shape,95); X1_05 = prctile(shape,5);
        X2_95 = prctile(scale_t,95); X2_05 = prctile(scale_t,5);

        for i = 1:length(Return_Period)
            RLtemp = gaminv(1 - 1/Return_Period(i), exp(X1), exp(X2));
            if no == 1
                CI_05 = gaminv(1 - 1/Return_Period(i), exp(X1_05), exp(X2_05));
                CI_95 = gaminv(1 - 1/Return_Period(i), exp(X1_95), exp(X2_95));
            end
            p(i) = plot(X_fit, RLtemp, '-.', 'Color', RP_COLORS(i,:), 'LineWidth', 1.4);
            RL_all{i} = RLtemp;
        end

    case 'lognormal'
        shape = params{1}; scale_t = params{2};
        X1 = mean(shape); X2 = mean(scale_t);
        X1_95 = prctile(shape,95); X1_05 = prctile(shape,5);
        X2_95 = prctile(scale_t,95); X2_05 = prctile(scale_t,5);

        for i = 1:length(Return_Period)
            RLtemp = logninv(1 - 1/Return_Period(i), X1, exp(X2));
            if no == 1
                CI_05 = logninv(1 - 1/Return_Period(i), X1_05, exp(X2_05));
                CI_95 = logninv(1 - 1/Return_Period(i), X1_95, exp(X2_95));
            end
            p(i) = plot(X_fit, RLtemp, '-.', 'Color', RP_COLORS(i,:), 'LineWidth', 1.4);
            RL_all{i} = RLtemp;
        end

    case 'weibull'
        shape = params{1}; scale_t = params{2};
        X1 = mean(shape); X2 = mean(scale_t);
        X1_95 = prctile(shape,95); X1_05 = prctile(shape,5);
        X2_95 = prctile(scale_t,95); X2_05 = prctile(scale_t,5);

        for i = 1:length(Return_Period)
            RLtemp = wblinv(1 - 1/Return_Period(i), exp(X1), exp(X2));
            if no == 1
                CI_05 = wblinv(1 - 1/Return_Period(i), exp(X1_05), exp(X2_05));
                CI_95 = wblinv(1 - 1/Return_Period(i), exp(X1_95), exp(X2_95));
            end
            p(i) = plot(X_fit, RLtemp, '-.', 'Color', RP_COLORS(i,:), 'LineWidth', 1.4);
            RL_all{i} = RLtemp;
        end

    case 'ev'
        shape = params{1}; scale_t = params{2};
        X1 = mean(shape); X2 = mean(scale_t);
        X1_95 = prctile(shape,95); X1_05 = prctile(shape,5);
        X2_95 = prctile(scale_t,95); X2_05 = prctile(scale_t,5);

        for i = 1:length(Return_Period)
            RLtemp = evinv(1 - 1/Return_Period(i), X1, exp(X2));
            if no == 1
                CI_05 = evinv(1 - 1/Return_Period(i), X1_05, exp(X2_05));
                CI_95 = evinv(1 - 1/Return_Period(i), X1_95, exp(X2_95));
            end
            p(i) = plot(X_fit, RLtemp, '-.', 'Color', RP_COLORS(i,:), 'LineWidth', 1.4);
            RL_all{i} = RLtemp;
        end

    case 'normal'
        shape = params{1}; scale_t = params{2};
        X1 = mean(shape); X2 = mean(scale_t);
        X1_95 = prctile(shape,95); X1_05 = prctile(shape,5);
        X2_95 = prctile(scale_t,95); X2_05 = prctile(scale_t,5);

        for i = 1:length(Return_Period)
            RLtemp = norminv(1 - 1/Return_Period(i), X1, exp(X2));
            if no == 1
                CI_05 = norminv(1 - 1/Return_Period(i), X1_05, exp(X2_05));
                CI_95 = norminv(1 - 1/Return_Period(i), X1_95, exp(X2_95));
            end
            p(i) = plot(X_fit, RLtemp, '-.', 'Color', RP_COLORS(i,:), 'LineWidth', 1.4);
            RL_all{i} = RLtemp;
        end


    otherwise
        error('Unknown best_fit: %s. Use ''gev'', ''gp'', ''gamma'' or ''lognormal''.', best_fit);
end

% ---- Plot AMS (solid data line) ----
pData = plot(YEAR_Range, AMS, '-', 'Color', dataColor, 'LineWidth', 1.8);



% ---- Titles & Labels ----
if no1 == 1
    title('Non-Stationary', 'FontWeight', 'bold');
elseif no1 == 2
    title('Stationary', 'FontWeight', 'bold');
end

% Only add ylabel for left column tiles (1,3,5)
if any(no1 == [3,4,6])
    ylabel(Label, 'FontName', 'Times New Roman', 'FontSize', 12);
end

% xlabel for all bottom row tiles if needed
if no3 == 1
    xlabel('Year', 'FontName', 'Times New Roman', 'FontSize', 12);
end



% ---- Legend ----
% if no2 == 1
%     legendHandles = [pData, p(1:6)];
%     legendLabels  = {'Data', '2-Year RP', '5-Year RP', '10-Year RP', '25-Year RP', '50-Year RP', '100-Year RP'};
%     hLegend = legend(legendHandles, legendLabels, 'Location', 'southoutside', 'Orientation', 'horizontal');
%     set(hLegend, 'FontSize', 9, 'FontName', 'Times New Roman', 'Box', 'off');
% end

% ---- Axes ----
xlim([min(X_fit) max(X_fit)]);
maxRL = max(cellfun(@(c) max(c(:)), RL_all));
minAMS = min(AMS(:));
ylim([minAMS * 0.9, maxRL * 1.1]);
grid on; box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);

% ---- Outputs ----
[RL, RL1, RL2, RL3, RL4, RL5] = deal(RL_all{1}, RL_all{2}, RL_all{3}, RL_all{4}, RL_all{5}, RL_all{6});

end

