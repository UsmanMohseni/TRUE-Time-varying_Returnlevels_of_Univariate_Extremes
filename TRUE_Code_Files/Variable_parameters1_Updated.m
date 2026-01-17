%% Function Variable Parameters
function [int_var]= Variable_parameters1_Updated(int_var,Model_type,change_point,AMS,par,Covariate,best_fit)
Years=size(AMS,1);
Covariate1=Covariate(11:Years-9,:);

[p,~]=size(Model_type);

        switch best_fit
            case 'gev'
                for i=1:1:p
                    if(strcmp(Model_type{i,1},'Stationary'))
                        p=gevfit(AMS);
                        
                        if(i==2)
                            int_var(1,2*i-1)=log(p(:,i));
                        else
                            int_var(1,2*i-1)=p(:,i);
                        end
                   elseif(strcmp(Model_type{i,1},'Non-Stationary Trend'))
                         Model=fitlm(Covariate1,par(i,:)');
                         coefficients = Model.Coefficients.Estimate;
                         B=[coefficients(2:end,1);coefficients(1,1)];
                         coefficients=B;
                         % int_var(:,2*i-1)=flipud(coefficients);
                         int_var(:,2*i-1)=coefficients;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=gevfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                         
                        p1=gevfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end,1),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=gevfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=gevfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    end
                end

            case 'gp'
                for i=1:1:p
                    if(strcmp(Model_type{i,1},'Stationary'))
                        p=gpfit(AMS);
                        if(i==2)
                            int_var(1,2*i-1)=log(p(:,i));
                        else
                            int_var(1,2*i-1)=p(:,i);
                        end
                      
                   elseif(strcmp(Model_type{i,1},'Non-Stationary Trend'))
                        Model=fitlm(Covariate1,par(i,:)');
                         coefficients = Model.Coefficients.Estimate;
                         B=[coefficients(2:end,1);coefficients(1,1)];
                         coefficients=B;
                         % int_var(:,2*i-1)=flipud(coefficients);
                         int_var(:,2*i-1)=coefficients;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=gpfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                         
                        p1=gpfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end,1),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=gpfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=gpfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    end
                end

            case 'lognormal'
                   for i=1:1:p
                    if(strcmp(Model_type{i,1},'Stationary'))
                        p=lognfit(AMS);
                        if(i==2)
                            int_var(1,2*i-1)=log(p(:,i));
                        else
                            int_var(1,2*i-1)=p(:,i);
                        end
                   elseif(strcmp(Model_type{i,1},'Non-Stationary Trend'))
                        Model=fitlm(Covariate1,par(i,:)');
                         coefficients = Model.Coefficients.Estimate;
                         B=[coefficients(2:end,1);coefficients(1,1)];
                         coefficients=B;
                         % int_var(:,2*i-1)=flipud(coefficients);
                         int_var(:,2*i-1)=coefficients;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=lognfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                         
                        p1=lognfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end,1),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=lognfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=lognfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    end
                   end

            case 'gamma'
              for i=1:1:p
                    if(strcmp(Model_type{i,1},'Stationary'))
                        p=gamfit(AMS);
                        int_var(1,2*i-1)=log(p(:,i));
                   elseif(strcmp(Model_type{i,1},'Non-Stationary Trend'))
                       Model=fitlm(Covariate1,par(i,:)');
                         coefficients = Model.Coefficients.Estimate;
                         B=[coefficients(2:end,1);coefficients(1,1)];
                         coefficients=B;
                         % int_var(:,2*i-1)=flipud(coefficients);
                         int_var(:,2*i-1)=coefficients;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=gamfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                         
                        p1=gamfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end,1),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=gamfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=gamfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    end
              end

              case 'normal'
                   for i=1:1:p
                    if(strcmp(Model_type{i,1},'Stationary'))
                        p=normfit(AMS);
                        if(i==2)
                            int_var(1,2*i-1)=log(p(:,i));
                        else
                            int_var(1,2*i-1)=p(:,i);
                        end
                   elseif(strcmp(Model_type{i,1},'Non-Stationary Trend'))
                        Model=fitlm(Covariate1,par(i,:)');
                         coefficients = Model.Coefficients.Estimate;
                         B=[coefficients(2:end,1);coefficients(1,1)];
                         coefficients=B;
                         % int_var(:,2*i-1)=flipud(coefficients);
                         int_var(:,2*i-1)=coefficients;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=normfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                         
                        p1=normfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end,1),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=normfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=normfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    end
                   end

            case 'ev'
                   for i=1:1:p
                    if(strcmp(Model_type{i,1},'Stationary'))
                        p=evfit(AMS);
                        if(i==2)
                            int_var(1,2*i-1)=log(p(:,i));
                        else
                            int_var(1,2*i-1)=p(:,i);
                        end
                   elseif(strcmp(Model_type{i,1},'Non-Stationary Trend'))
                        Model=fitlm(Covariate1,par(i,:)');
                         coefficients = Model.Coefficients.Estimate;
                         B=[coefficients(2:end,1);coefficients(1,1)];
                         coefficients=B;
                         % int_var(:,2*i-1)=flipud(coefficients);
                         int_var(:,2*i-1)=coefficients;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=evfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                         
                        p1=evfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end,1),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=evfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=evfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    end
                   end

            case 'Weibull'
              for i=1:1:p
                    if(strcmp(Model_type{i,1},'Stationary'))
                        p=wblfit(AMS);
                        int_var(1,2*i-1)=log(p(:,i));
                   elseif(strcmp(Model_type{i,1},'Non-Stationary Trend'))
                       Model=fitlm(Covariate1,par(i,:)');
                         coefficients = Model.Coefficients.Estimate;
                         B=[coefficients(2:end,1);coefficients(1,1)];
                         coefficients=B;
                         % int_var(:,2*i-1)=flipud(coefficients);
                         int_var(:,2*i-1)=coefficients;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=polyfit([11:10+change_point],par(i,1:change_point),1);
                        int_var(:,2*i-1)=p1;
                        p2=wblfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step Trend'))
                         
                        p1=wblfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=polyfit([11:10+change_point],par(i,change_point+1:end,1),1);
                        int_var(:,2*i)=p2;
                    elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
                        p1=wblfit(AMS(1:change_point,1));
                        int_var(1,2*i-1)=p1(:,i);
                        p2=wblfit(AMS(change_point+1:end,1));
                        int_var(1,2*i)=p2(:,i);
                    end
              end

        case 'logistic'
            for i = 1:1:p
        
                if strcmp(Model_type{i,1},'Stationary')
                    pd = fitdist(AMS(:,1),'Logistic');
        
                    if i == 2
                        int_var(1,2*i-1) = log(pd.sigma);   % scale
                    else
                        int_var(1,2*i-1) = pd.mu;           % location
                    end
        
                elseif strcmp(Model_type{i,1},'Non-Stationary Trend')
                    Model = fitlm(Covariate1, par(i,:)');
                    coefficients = Model.Coefficients.Estimate;
                    B = [coefficients(2:end,1); coefficients(1,1)];
                    int_var(:,2*i-1) = B;
        
                elseif strcmp(Model_type{i,1},'Step Trend') && strcmp(Model_type{i,2},'Step Trend')
                    p1 = polyfit([11:10+change_point], par(i,1:change_point), 1);
                    int_var(:,2*i-1) = p1;
        
                    p2 = polyfit([11:10+change_point], par(i,change_point+1:end), 1);
                    int_var(:,2*i) = p2;
        
                elseif strcmp(Model_type{i,1},'Step Trend') && strcmp(Model_type{i,2},'Step No-Trend')
                    p1 = polyfit([11:10+change_point], par(i,1:change_point), 1);
                    int_var(:,2*i-1) = p1;
        
                    pd2 = fitdist(AMS(change_point+1:end,1),'Logistic');
                    if i == 2
                        int_var(1,2*i) = log(pd2.sigma);
                    else
                        int_var(1,2*i) = pd2.mu;
                    end
        
                elseif strcmp(Model_type{i,1},'Step No-Trend') && strcmp(Model_type{i,2},'Step Trend')
                    pd1 = fitdist(AMS(1:change_point,1),'Logistic');
                    if i == 2
                        int_var(1,2*i-1) = log(pd1.sigma);
                    else
                        int_var(1,2*i-1) = pd1.mu;
                    end
        
                    p2 = polyfit([11:10+change_point], par(i,change_point+1:end,1), 1);
                    int_var(:,2*i) = p2;
        
                elseif strcmp(Model_type{i,1},'Step No-Trend') && strcmp(Model_type{i,2},'Step No-Trend')
                    pd1 = fitdist(AMS(1:change_point,1),'Logistic');
                    if i == 2
                        int_var(1,2*i-1) = log(pd1.sigma);
                    else
                        int_var(1,2*i-1) = pd1.mu;
                    end
        
                    pd2 = fitdist(AMS(change_point+1:end,1),'Logistic');
                    if i == 2
                        int_var(1,2*i) = log(pd2.sigma);
                    else
                        int_var(1,2*i) = pd2.mu;
                    end
                    end
                end

        end
end
