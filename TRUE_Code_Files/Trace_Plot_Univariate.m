function Trace_Plot_Univariate(para, Model_type, int_var, best_fit, No)
    [m, n] = size(para);
    z = 1;  % Subplot index
    y = 1;  % Parameter index
    if (n > 3)
        X = 16.88;
    else
        X = 8.5;
    end

    % Create the figure
    figure(No);
    set(gcf, 'Color', 'w', 'Units', 'centimeters', 'Position', [9,1.05, 19.55, X]); % Set size to 8x10 inches

    % Colors
    traceColor = [0, 0.6, 0.5];
    histColor = [1, 0.6, 0.2];
    fitColor = [0.5, 0, 0.5];
    gridColor = [0.8, 0.8, 0.8];

    % Parameter names
    switch best_fit
        case 'gev'
            param_names = {'\gamma', '\sigma', '\mu'};
        case 'gp'
            param_names = {'\gamma', '\sigma'};
        case 'gamma'
            param_names = {'\alpha', '\beta'};
        case 'lognormal'
            param_names = {'\mu', '\sigma'};
        case 'normal'
            param_names = {'\mu', '\sigma'};
        case 'ev'
            param_names = {'\mu', '\sigma'};
        case 'weibull'
            param_names = {'\alpha', '\beta'};
        case 'logistic'
            param_names = {'\mu', '\sigma'};
    end

    total_rows = n; % rows of subplot
    middle_hist_idx = ceil(total_rows / 2); % middle histogram row

    % Loop
    for i = 1:size(Model_type, 1)
        model_condition = Model_type(i, 1);

        if strcmp(model_condition, 'Stationary')
            % Trace plot
            subplot(n, 2, z);
            plot(para(:, y), 'Color', traceColor, 'LineWidth', 1.5);
            axis tight; grid on;
            set(gca, 'GridColor', gridColor);

            % Y-axis label for all trace plots
            ylabel(param_names{i}, 'FontWeight', 'bold', 'FontSize', 12);

            % X-axis label only for bottom trace plot
            if z == (total_rows * 2 - 1)
                xlabel('Iteration', 'FontWeight', 'bold', 'FontSize', 12);
            else
                xlabel('');
            end

            % Histogram plot
            z = z + 1;
            subplot(n, 2, z);
            hold on;
            pd = fitdist(para(:, y), 'Normal');
            x = linspace(min(para(:, y)), max(para(:, y)), 100);
            y_fit = pdf(pd, x);
            histogram(para(:, y), 100, 'Normalization', 'pdf', ...
                'FaceColor', histColor, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            plot(x, y_fit, 'Color', fitColor, 'LineWidth', 2);

            % Y-axis label only for middle histogram plot
            if (z/2) == middle_hist_idx
                ylabel('Probability Density', 'FontWeight', 'bold', 'FontSize', 12);
            else
                ylabel('');
            end

            % X-axis label for all histograms
            xlabel(param_names{i}, 'FontWeight', 'bold', 'FontSize', 12);

            z = z + 1;
            y = y + 1;

        elseif strcmp(model_condition, 'Non-Stationary Trend')
            for j = 1:size(int_var, 1)
                param_label = [param_names{i}, '_', num2str(j-1)];

                % Trace plot
                subplot(n, 2, z);
                plot(para(:, y), 'Color', traceColor, 'LineWidth', 1.5);
                axis tight; grid on;
                set(gca, 'GridColor', gridColor);

                % Y-axis for all trace plots
                ylabel(param_label, 'FontWeight', 'bold', 'FontSize', 12);

                % X-axis only for bottom trace plot
                if z == (total_rows * 2 - 1)
                    xlabel('Iteration', 'FontWeight', 'bold', 'FontSize', 12);
                else
                    xlabel('');
                end

                % Histogram plot
                z = z + 1;
                subplot(n, 2, z);
                hold on;
                pd = fitdist(para(:, y), 'Normal');
                x = linspace(min(para(:, y)), max(para(:, y)), 100);
                y_fit = pdf(pd, x);
                histogram(para(:, y), 100, 'Normalization', 'pdf', ...
                    'FaceColor', histColor, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                plot(x, y_fit, 'Color', fitColor, 'LineWidth', 2);

                % Y-axis only for middle histogram plot
                if (z/2) == middle_hist_idx
                    ylabel('Probability Density', 'FontWeight', 'bold', 'FontSize', 12);
                else
                    ylabel('');
                end

                % X-axis for all histograms
                xlabel(param_label, 'FontWeight', 'bold', 'FontSize', 12);

                z = z + 1;
                y = y + 1;
            end
        end
    end

    sgtitle('Parameter Trace and Histogram Plots', 'FontWeight', 'bold', 'FontSize', 16);
end

