function Univariate_Parameter_Trend_Plot(int_MW_par1, dist_type, years, MW, No)
    % Adjust years to match moving window
    years = years(MW/2+1:end-MW/2);

    % Get the size for the first distribution only
    [~, NN] = size(int_MW_par1(:, :, 1));

    % Create figure
    figure(No);
    set(gcf, 'Color', 'w', 'Units', 'centimeters', 'Position', [12.3, 4.9, 12.75, 6.3]); % Set size to 8x10 inches

    % Colors
    vibrantColor = [0.1, 0.62, 0.94];  
    trendColor   = [0.85, 0.33, 0.1];  

    % Use only the first distribution
    j = 1; 
    switch lower(dist_type{j})
            case 'gev'
                y_label = {'\sigma', '\mu'};
                param_idx = [2, 3]; % Skip shape (k)
                scale_param_row = 2; % σ is 2nd row
            case 'gp'
                y_label = {'\sigma', '\gamma'};
                param_idx = [1, 2];
                scale_param_row = 1; % σ is 1st row
            case 'gamma'
                y_label = {'\alpha', '\beta'};
                param_idx = [1, 2];
                scale_param_row = 2; % β is 2nd row
            case 'lognormal'
                y_label = {'\mu', '\sigma'};
                param_idx = [1, 2];
                scale_param_row = 2; % σ is 2nd row
            case 'normal'
                y_label = {'\mu', '\sigma'};
                param_idx = [1, 2];
                scale_param_row = 2; % σ is 2nd row
            case 'ev'
                y_label = {'\mu', '\sigma'};
                param_idx = [1, 2];
                scale_param_row = 2; % σ is 2nd row
            case 'logistic'
                y_label = {'\mu', '\sigma'};
                param_idx = [1, 2];
                scale_param_row = 2; % σ is 2nd row
            case 'weibull'
                y_label = {'\alpha', '\beta'};
                param_idx = [1, 2];
                scale_param_row = 2; % β is 2nd row
        end

    % Extract only the first distribution parameters
    int_MW_par = int_MW_par1(:, :, j);

    % Fix negative scale parameter
        if any(int_MW_par(scale_param_row, :) <= 0)
            warning('Negative or zero scale parameter detected for %s. Transforming...', dist_type{j});
            int_MW_par(scale_param_row, :) = abs(int_MW_par(scale_param_row, :) + eps);
        end

    % 1×2 tiled layout
    tiledlayout(1, 2, "TileSpacing", "loose", "Padding", "compact");

    % Loop over the two parameters
    for i = 1:2
        nexttile;

        % Plot data
        plot(int_MW_par(i, :), 'Color', vibrantColor, 'LineWidth', 1.5);
        hold on;

        % Fit trend
        Coeff = polyfit(1:NN, int_MW_par(i, :), 1);
        Trend = polyval(Coeff, 1:NN);
        plot(1:NN, Trend, '-.', 'Color', trendColor, 'LineWidth', 2);

        % Calculate slope
        S = Coeff(1);

        % Add slope text in top-left corner
        text(0.05, 0.9, sprintf('S = %.4f', S), ...
            'Units', 'normalized', 'FontName', 'Times New Roman', ...
            'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');

        % Axis limits
        xlim([-5, NN + 5]);
        ylim([min(int_MW_par(i, :)) - 0.5, max(int_MW_par(i, :)) + 0.5]);

        % X ticks & labels
        xticks(1:10:NN);
        xticklabels(years(1:10:end));
        xtickangle(45);

        % Labels
        xlabel('Year', 'FontName', 'Times New Roman', 'FontSize', 8, 'FontWeight', 'bold');
        ylabel(y_label{i}, 'FontName', 'Times New Roman', 'FontSize', 8, 'FontWeight', 'bold');

        % Style
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'Box', 'on', ...
            'TickDir', 'out', 'LineWidth', 1);

        hold off;
    end
end

