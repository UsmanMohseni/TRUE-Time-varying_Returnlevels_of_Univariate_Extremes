function [RL, RL1, RL2]=Return_level_plot(params,AMS,YEAR_Range,best_fit,Label,no,no1,no2,no3,no4)
% no= To show CI or not
% no1= anything related to first tile
% no2=anything related to second tile
% no3 =anything related to third tile
COLOR=[0.5 0 0;0 0.5 0;0.9290 0.6940 0.1250];
band_color=[1 0.7 0.7; 0.7 1 0.7;1 1 0.7];
Return_Period=[10;40;100];
% Set figure background to white
set(gcf, 'Color', 'white'); 

% After your plot commands, set axes background to white
set(gca, 'Color', 'white');

 nexttile

switch best_fit
    case 'gev'
        shape=params{1};
        scale_t=params{2};
        mu_t=params{3};

        X1=mean(shape);
        X1_95=(prctile(shape,95));
        X1_05=(prctile(shape,5));
        X2=mean(scale_t);
        X2_95=(prctile(scale_t,95));
        X2_05=(prctile(scale_t,5));
        X3=mean(mu_t);
        X3_95=(prctile(mu_t,95));
        X3_05=(prctile(mu_t,5));
        X_fit=YEAR_Range;
        
        RL=(gevinv(1-1/Return_Period(1,1),X1,exp(X2),X3));
        CI_05=(gevinv(1-1/Return_Period(1,1),X1_05,exp(X2_05),X3_05));

        CI_95=(gevinv(1-1/Return_Period(1,1),X1_95,exp(X2_95),X3_95));
if(no==1)
        
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(1,:),'edgecolor','none')
             hold on
             end 
        plot(X_fit,RL,'Color',COLOR(1,:),'LineWidth',2)
            hold on
   
        
         RL1=gevinv(1-1/Return_Period(2,1),X1,exp(X2),X3);
         CI_05=(gevinv(1-1/Return_Period(2,1),X1_05,exp(X2_05),X3_05));
  
        CI_95=(gevinv(1-1/Return_Period(2,1),X1_95,exp(X2_95),X3_95));

if(no==1)        
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(2,:),'edgecolor','none')
             hold on
end     
            plot(X_fit,RL1,'Color',COLOR(2,:),'LineWidth',2)
            hold on

          RL2=gevinv(1-1/Return_Period(3,1),X1,exp(X2),X3);
           CI_05=(gevinv(1-1/Return_Period(3,1),X1_05,exp(X2_05),X3_05));

        CI_95=(gevinv(1-1/Return_Period(3,1),X1_95,exp(X2_95),X3_95));
if(no==1)   
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(3,:),'edgecolor','none')
             hold on
end
            plot(X_fit,RL2,'color',COLOR(3,:),'LineWidth',2)
            hold on

        
            plot(YEAR_Range,AMS,'blue','LineWidth',3)
        hold on
        %ylim=[0, max(RL2)*1.5];
    case 'gp'
        shape=params{1};
        scale_t=params{2};
       

        X1=mean(shape);
        X1_95=(prctile(shape,95));
        X1_05=(prctile(shape,5));
        X2=mean(scale_t);
        X2_95=(prctile(scale_t,95));
        X2_05=(prctile(scale_t,5));
       
        X_fit=YEAR_Range;
        
        RL=(gpinv(1-1/Return_Period(1,1),X1,exp(X2)));
        CI_05=(gpinv(1-1/Return_Period(1,1),X1_05,exp(X2_05)));

        CI_95=(gpinv(1-1/Return_Period(1,1),X1_95,exp(X2_95)));

        if(no==1)
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(1,:),'edgecolor','none')
             hold on
              end 
        plot(X_fit,RL,'Color',COLOR(1,:),'LineWidth',2)
            hold on
         
        
         RL1=gpinv(1-1/Return_Period(2,1),X1,exp(X2));
         CI_05=(gpinv(1-1/Return_Period(2,1),X1_05,exp(X2_05)));

        CI_95=(gpinv(1-1/Return_Period(2,1),X1_95,exp(X2_95)));

        if(no==1)
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(2,:),'edgecolor','none')
             hold on
        end
            plot(X_fit,RL1,'Color',COLOR(2,:),'LineWidth',2)
            hold on
       
          RL2=gpinv(1-1/Return_Period(3,1),X1,exp(X2));
           CI_05=(gpinv(1-1/Return_Period(3,1),X1_05,exp(X2_05)));
        CI_05=(CI_05(1:size(X1,2)));
        CI_95=(gpinv(1-1/Return_Period(3,1),X1_95,exp(X2_95)));
        CI_95=(CI_95(1:size(X1,2)));
         RL2=(RL2(1:size(X1,2)));
        if(no==1)
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(3,:),'edgecolor','none')
             hold on
        end
            plot(X_fit,RL2,'color',COLOR(3,:),'LineWidth',2)
            hold on
        
            plot(YEAR_Range,AMS,'blue','LineWidth',3)
        hold on

        % ylim=[0, max(RL2)*1.5];
    case 'gamma'
        shape=params{1};
        scale_t=params{2};
       

        X1=mean(shape);
        X1_95=(prctile(shape,95));
        X1_05=(prctile(shape,5));
        X2=mean(scale_t);
        X2_95=(prctile(scale_t,95));
        X2_05=(prctile(scale_t,5));
       
        X_fit=YEAR_Range;
        
        RL=(gaminv(1-1/Return_Period(1,1),exp(X1),exp(X2)));
        CI_05=(gaminv(1-1/Return_Period(1,1),exp(X1_05),exp(X2_05)));

        CI_95=(gaminv(1-1/Return_Period(1,1),exp(X1_95),exp(X2_95)));

        if(no==1)
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(1,:),'edgecolor','none')
             hold on
        end
        plot(X_fit,RL,'Color',COLOR(1,:),'LineWidth',2)
            hold on
            
     
         RL1=gaminv(1-1/Return_Period(2,1),exp(X1),exp(X2));
         CI_05=(gaminv(1-1/Return_Period(2,1),exp(X1_05),exp(X2_05)));
      
        CI_95=(gaminv(1-1/Return_Period(2,1),exp(X1_95),exp(X2_95)));
       
          if(no==1)
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(2,:),'edgecolor','none')
             hold on
          end
            plot(X_fit,RL1,'Color',COLOR(2,:),'LineWidth',2)
            hold on
          RL2=gaminv(1-1/Return_Period(3,1),exp(X1),exp(X2));
           CI_05=(gaminv(1-1/Return_Period(3,1),exp(X1_05),exp(X2_05)));
   
        CI_95=(gaminv(1-1/Return_Period(3,1),exp(X1_95),exp(X2_95)));
      
          if(no==1)
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(3,:),'edgecolor','none')
             hold on
          end
            plot(X_fit,RL2,'color',COLOR(3,:),'LineWidth',2)
            hold on
        
            plot(YEAR_Range,AMS,'blue','LineWidth',3)
        hold on
        % ylim[0, max(RL2)*1.5];
    case 'lognormal'
        shape=params{1};
        scale_t=params{2};
       

        X1=mean(shape);
        X1_95=(prctile(shape,95));
        X1_05=(prctile(shape,5));
        X2=mean(scale_t);
        X2_95=(prctile(scale_t,95));
        X2_05=(prctile(scale_t,5));
       
        X_fit=YEAR_Range;
        
        RL=(logninv(1-1/Return_Period(1,1),X1,exp(X2)));
        CI_05=(logninv(1-1/Return_Period(1,1),X1_05,exp(X2_05)));
        CI_05=(CI_05(1:size(X1,2)));
        CI_95=(logninv(1-1/Return_Period(1,1),X1_95,exp(X2_95)));
        CI_95=(CI_95(1:size(X1,2)));
        RL=(RL(1:size(X1,2))) ;
          if(no==1)
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(1,:),'edgecolor','none')
             hold on
          end
        plot(X_fit,RL,'Color',COLOR(1,:),'LineWidth',2)
            hold on
            
        
         RL1=logninv(1-1/Return_Period(2,1),X1,exp(X2));
         CI_05=(logninv(1-1/Return_Period(2,1),X1_05,exp(X2_05)));
        CI_05=(CI_05(1:size(X1,2)));
        CI_95=(logninv(1-1/Return_Period(2,1),X1_95,exp(X2_95)));
        CI_95=(CI_95(1:size(X1,2)));
         RL1=(RL1(1:size(X1,2)));
          if(no==1)
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(2,:),'edgecolor','none')
             hold on
          end
            plot(X_fit,RL1,'Color',COLOR(2,:),'LineWidth',2)
            hold on
          RL2=logninv(1-1/Return_Period(3,1),X1,exp(X2));
           CI_05=(logninv(1-1/Return_Period(3,1),X1_05,exp(X2_05)));
        CI_05=(CI_05(1:size(X1,2)));
        CI_95=(logninv(1-1/Return_Period(3,1),X1_95,exp(X2_95)));
        CI_95=(CI_95(1:size(X1,2)));
         RL2=(RL2(1:size(X1,2)));
          if(no==1)
        %fill([X_fit',fliplr(X_fit')],[CI_05,fliplr(CI_95)],band_color(3,:),'edgecolor','none')
             hold on
          end
            plot(X_fit,RL2,'color',COLOR(3,:),'LineWidth',2)
            hold on
        
            plot(YEAR_Range,AMS,'blue','LineWidth',3)
        hold on
        % ylim=[0, max(RL2)*1.5];
end


if(no1==1)
% After your plotting commands
 title("Non-Stationary");
elseif(no1==2)

 title("Stationary");
end
if(no2==1)
    
hLegend=legend('10Yr CI','10Yr RL','40Yr CI','40Yr RL','100Yr CI','100Yr RL','AMS','Location','southoutside','Orientation','horizontal')
 set(hLegend,'Position', [0.038098965046617,0.020061728395062,0.915364569906767,0.024305554969167]); % Adjust values as needed
end
if(no3==1)
xlabel('Year')
end
if(no4==1)
ylabel(Label)
end
xlim([min(X_fit) max(X_fit)])
ylim([0 max(RL2)*1.2])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
end
