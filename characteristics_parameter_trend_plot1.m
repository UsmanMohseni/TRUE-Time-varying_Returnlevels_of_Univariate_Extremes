function characteristics_parameter_trend_plot1(int_MW_par1, dist_type,No)
    % Define labels, legends, and titles based on distribution type
    switch lower(dist_type)
        case 'gev'
            LEG = {'Scale Parameter', 'Trend'; 'Location Parameter', 'Trend'};
            TITLE = {'Trend in Scale Parameter', 'Trend in Location Parameter'};
            y_label = { '\sigma', '\mu'};
           

        case 'gp'
            LEG = {'Scale Parameter', 'Trend'; 'Shape Parameter', 'Trend'};
            TITLE = {'Trend in Scale Parameter', 'Trend in Shape Parameter'};
            y_label = {'\sigma', '\gamma'};
            

        case 'gamma'
             LEG = {'Alpha Parameter', 'Trend'; 'Beta Parameter', 'Trend'};
            TITLE = {'Trend in Alpha Parameter', 'Trend in Beta Parameter'};
            y_label = {'\alpha', '\beta'};
            
        case 'lognormal'
              LEG = {'Mean', 'Trend'; 'Standard deviation', 'Trend'};
            TITLE = {'Trend in Alpha Parameter', 'Trend in Beta Parameter'};
            y_label = {'\alpha', '\beta'};
           
    end

    % Get the size of the input data
    [~, NN] = size(int_MW_par1(:, :, 1));
    figure(No)
    
    % Plot index to arrange subplots (6-panel layout)
    plotindex = [1 4 2 5 3 6];
    ind = 1;

    % Loop over parameters based on the distribution type
    for j = 1:3
          switch lower(dist_type{j,1})
        case 'gev'
            LEG = {'Scale Parameter', 'Trend'; 'Location Parameter', 'Trend'};
            TITLE = {'Trend in Scale Parameter', 'Trend in Location Parameter'};
            y_label = { '\sigma', '\mu'};
           

        case 'gp'
            LEG = {'Scale Parameter', 'Trend'; 'Shape Parameter', 'Trend'};
            TITLE = {'Trend in Scale Parameter', 'Trend in Shape Parameter'};
            y_label = {'\sigma', '\gamma'};
            

        case 'gamma'
             LEG = {'Alpha Parameter', 'Trend'; 'Beta Parameter', 'Trend'};
            TITLE = {'Trend in Alpha Parameter', 'Trend in Beta Parameter'};
            y_label = {'\alpha', '\beta'};
            
        case 'lognormal'
              LEG = {'Mean', 'Trend'; 'Standard deviation', 'Trend'};
            TITLE = {'Trend in Alpha Parameter', 'Trend in Beta Parameter'};
            y_label = {'\alpha', '\beta'};
           
    end

        int_MW_par = int_MW_par1(:, :, j);  % Extract parameter-specific data
        
        for i = 1:2  % Loop over scale and shape (or relevant) parameters
            subplot(2, 3, plotindex(ind))
            plot(int_MW_par(i, :), 'color', 'black')
            
            % Trend line using polyfit
            Coeff = polyfit(1:NN, int_MW_par(i, :), 1);
            Trend = polyval(Coeff, 1:NN);
            xlim([1, size(int_MW_par, 2)])
            ylim([min(int_MW_par(i, :)), max(int_MW_par(i, :))])
            hold on
            plot(1:NN, Trend, 'r-.', 'LineWidth', 2);
            hold on

            % Set labels and titles dynamically
            xlabel('Moving Window');
            ylabel(y_label{i });
            legend(LEG(i , :), 'Location', 'southoutside', 'Orientation', 'horizontal')
            title(TITLE{i})

            ind = ind + 1;
        end
    end
end
