clc
clear
close all
load D:\MTech\Thesis_RESULTS\India\Narmada_Basin\MATLAB_FILES\Comparison_plot_files\EMD_AP.mat
RL_STA=RL;
RL1_STA=RL1;
RL2_STA=RL2;
load D:\MTech\Thesis_RESULTS\India\Narmada_Basin\MATLAB_FILES\Comparison_plot_files\NON_EMD_AP.mat
RL_NON_STA=RL;
RL1_NON_STA=RL1;
RL2_NON_STA=RL2;
Year=years;

figure(10)
plot(years,AMS,'-','color','black','LineWidth',2)
hold on
plot(years,RL_NON_STA','-.','color','black','LineWidth',1)
hold on
plot(years,RL_STA','--','color','red','LineWidth',1)
hold on
% plot(Year,RL1_NON_STA','-','color','black','LineWidth',1)
% hold on
% plot(Year,RL1_STA','--','color','black','LineWidth',1)
% hold on
% plot(Year,RL2_NON_STA','-.','color','red','LineWidth',1)
% hold on
% plot(Year,RL2_STA','--','color','black','LineWidth',1)
% hold on
xlabel('Year')
ylabel('Discharge (m^3,s)')
title('Time Varying Return Level Comparision Between EEMD and Non-EEMD Covariate Models')
legend('Peak Flow  Series','Model with Original Covariates','Model with EEMD covarites','Orientation','horizontal','Location','southoutside')