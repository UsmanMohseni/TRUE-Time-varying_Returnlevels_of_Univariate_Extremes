function parameter_trend_plot(int_MW_par,No)
    LEG={'Shape Parameter','Trend';'Scale Parameter','Trend';'Location Parameter','Trend'};
    TITLE={'Trend in Shape Parameter','Trend in Scale Parameter','Trend in Location Parameter'}; 
    y_label={'\gamma','\sigma','\mu'};
    [~,NN]=size(int_MW_par);
 figure(No)
for i=2:1:3
   
    subplot(2,1,i-1)
    plot(int_MW_par(i,:),'color','black')
    Coeff= polyfit(1:NN,int_MW_par(i,:),1);
    Trend=polyval(Coeff,1:NN);
    xlim([1,size(int_MW_par,2)])
    ylim([min(int_MW_par(i,:)),max(int_MW_par(i,:))])
    hold on
    plot(1:NN,Trend,'r-.','LineWidth',2);
    hold on
    xlabel('Moving Window');
    ylabel(y_label{1,i});
    legend(LEG(i,:),'Location','southoutside','Orientation','horizontal')
    title(TITLE{1,i})
end
% set(figure(2), 'WindowState', 'maximized');
end

