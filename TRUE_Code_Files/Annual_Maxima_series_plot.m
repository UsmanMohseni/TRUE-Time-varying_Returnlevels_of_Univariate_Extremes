function Annual_Maxima_series_plot(YEAR_Range,AMS)
   figure(1)
    plot(YEAR_Range,AMS,'b','LineWidth',2)
    % Coeff= polyfit(YEAR_Range,AMS,1);
    % Trend=polyval(Coeff,YEAR_Range);
    hold on
%     plot(YEAR_Range(1,1):YEAR_Range(1,2), Trend, 'r-', 'LineWidth', 2);
%     hold on
    xlabel('Year');
    ylabel('Peak Flow (m^3/sec)');
%     legend('Data','Trend','location','southoutside','Orientation','horizontal')
    title('Peak Flow Time Series')
    set(figure(1), 'WindowState', 'maximized');
end