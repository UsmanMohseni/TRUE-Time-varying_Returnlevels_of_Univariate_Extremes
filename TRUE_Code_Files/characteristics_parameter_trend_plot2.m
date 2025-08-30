function characteristics_parameter_trend_plot2(int_MW_par1, dist_type,years,MW,No)
years=years(MW/2+1:end-MW/2);
    % Get the size of the input data
    [~, NN] = size(int_MW_par1(:, :, 1));
    
    % Create figure with specified size (in inches) for Word embedding
    figure(No);
    set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1, 1, 5, 3.5]);  
    % [left, bottom, width, height] in inches

    % Plot index for 6-panel layout
    plotindex = [1, 4, 2, 5, 3, 6]; 
    ind = 1;

    % Define consistent colors
    vibrantColor = [0.1, 0.62, 0.94];  % Bright blue for data lines
    trendColor = [0.85, 0.33, 0.1];    % Red-orange for trend lines
TITLE={'Peak Flow ';'';'Duration';'';'Volume';''};
x=1;
tiledlayout(2,3,"TileSpacing","loose","Padding","tight")
    % Loop over parameters based on the distribution type
    for j = 1:3
        switch lower(dist_type{j})
            case 'gev'
                LEG = {'Scale Parameter', 'Trend'; 'Location Parameter', 'Trend'};
                % TITLE = {'Trend in Scale Parameter', 'Trend in Location Parameter'};
                y_label = {'\sigma', '\mu'};

            case 'gp'
                LEG = {'Scale Parameter', 'Trend'; 'Shape Parameter', 'Trend'};
                % TITLE = {'Trend in Scale Parameter', 'Trend in Shape Parameter'};
                y_label = {'\sigma', '\gamma'};

            case 'gamma'
                LEG = {'Alpha Parameter', 'Trend'; 'Beta Parameter', 'Trend'};
                % TITLE = {'Trend in Alpha Parameter', 'Trend in Beta Parameter'};
                y_label = {'\alpha', '\beta'};

            case 'lognormal'
                LEG = {'Mean', 'Trend'; 'Standard Deviation', 'Trend'};
                % TITLE = {'Trend in Alpha Parameter', 'Trend in Beta Parameter'};
                y_label = {'\alpha', '\beta'};
        end

        % Extract parameter-specific data
        int_MW_par = int_MW_par1(:, :, j);

        % Loop over two parameters (e.g., Scale and Location)
        for i = 1:2
            % subplot(2, 3, plotindex(ind));
            nexttile(plotindex(ind));
            
            % Plot the data line
            plot(int_MW_par(i, :), 'Color', vibrantColor, 'LineWidth', 1.5);
            hold on;

            % Fit and plot the trend line
            Coeff = polyfit(1:NN, int_MW_par(i, :), 1);
            Trend = polyval(Coeff, 1:NN);
            plot(1:NN, Trend, '-.', 'Color', trendColor, 'LineWidth', 2);

            % Adjust axis limits for zoomed-out view
            xlim([-5, NN + 5]);  % Extra space on x-axis
            ylim([min(int_MW_par(i, :)) - 0.5, max(int_MW_par(i, :)) + 0.5]);  % Extra padding on y-axis
            
            % Set custom xticks and labels to display years
            xticks(1:10:NN);  % Show xticks at intervals of 5 units
            xticklabels(years(1:10:end));  % Display the corresponding years   
            xtickangle(45)
            % Set labels and legend
            if(ind==2||ind==4||ind==6)
            xlabel('Year', 'FontName', 'Times New Roman', 'FontSize', 8,'FontWeight','bold');
            end
            ylabel(y_label{i}, 'FontName', 'Times New Roman', 'FontSize', 8,'FontWeight','bold');
            % legend(LEG(i, :), 'Location', 'southoutside', 'Orientation', 'vertical');

            % Add title
            title(TITLE{x}, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
x=x+1;
            % Improve plot appearance
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'Box', 'on', ...
                'TickDir', 'out', 'LineWidth', 1);
            hold off;

            % Increment subplot index
            ind = ind + 1;
        end
    end

end
