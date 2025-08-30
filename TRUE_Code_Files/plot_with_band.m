function plot_with_band(shape,scale_t,mu_t,sd1,mu1,Years,AMS,no)

%%plot with band
MED(1,:)=median(shape);
MED(2,:)=median(scale_t);
MED(3,:)=median(mu_t);


% percentile95(1,:) =prctile(shape,99,1);
% percentile95(2,:)=prctile(scale_t,99,1);
% percentile95(3,:)=prctile(mu_t,99,1);
% 
% percentile5(1,:) =prctile(shape,1,1);
% percentile5(2,:)=prctile(scale_t,1,1);
% percentile5(3,:)=prctile(mu_t,1,1);
% 
% percentile25(1,:) =prctile(shape,25,1);
% percentile25(2,:)=prctile(scale_t,25,1);
% percentile25(3,:)=prctile(mu_t,25,1);
% 
% percentile75(1,:) =prctile(shape,75,1);
% percentile75(2,:)=prctile(scale_t,75,1);
% percentile75(3,:)=prctile(mu_t,75,1);

T=2;
Title=sprintf(' %d Year Return Level (Discharge)',T);
RLmed=gevinv(0.5,MED(1,:),MED(2,:),MED(3,:)).*sd1+mu1;

RL95=gevinv(0.9,MED(1,:),MED(2,:),MED(3,:)).*sd1+mu1;

RL5=gevinv(0.1,MED(1,:),MED(2,:),MED(3,:)).*sd1+mu1;

RL25=gevinv(0.25,MED(1,:),MED(2,:),MED(3,:)).*sd1+mu1;

RL75=gevinv(0.75,MED(1,:),MED(2,:),MED(3,:)).*sd1+mu1;

% Specify the color for the band
color = [0 0.4470 0.7410;0 1 1]; % Gray color
figure (no)


% Plot the color band
fill([Years', fliplr(Years')], [RL95, fliplr(RL5)], color(2,:),'EdgeColor','none');
hold on
fill([Years', fliplr(Years')], [RL75, fliplr(RL25)], color(1,:),'EdgeColor','none');
hold on
plot(Years,RLmed,'color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
scatter (Years,AMS,25,'black','filled')
hold on
xlabel('Years')
ylabel('Discharge(Cumec)')
legend('5%/95%','25%/75%','Median (50%)','Peak FLow Observations','Location','southoutside','Orientation','horizontal')
%title (Title)

end