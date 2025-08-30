
function [DIC1,DIC2,PARA_FINAL]= model_evaluator2(numVariables,allCombinations,Covariate1,AMS,best_fit)
warning('off', 'all');
com=1;
PARA_FINAL(1:25005,1:15,1:63)=nan;


[~,Model_type1,~]=initial_var(AMS,Covariate1,best_fit,20);

if(any(strcmp(Model_type1(:,1),'Non-Stationary Trend')))

    for variableIndex = 1:numVariables
        combinationsForVariable = allCombinations{variableIndex};
        
        % Loop through each combination
        for combinationIndex = 1:size(combinationsForVariable, 1)
            currentCombination = combinationsForVariable(combinationIndex, :);
            if(currentCombination(:,end)~=7)
                 tic;
            % Extract data for the current combination
            Covariate = Covariate1(:, currentCombination);
           
            numPredictors = size(Covariate, 2);
    VIF = zeros(numPredictors, 1);
    
    for i = 1:numPredictors
        predictorsWithoutI = [Covariate(:, 1:i-1), Covariate(:, i+1:end)];
        mdlWithoutI = fitlm(predictorsWithoutI, Covariate(:, i));
        R2 = mdlWithoutI.Rsquared.Ordinary;
        VIF(i) = 1 / (1 - R2);
    end
    index=0; % For checking whether we need PCA or not
    for i=1:1:length(VIF)
        if(VIF(i,1)>3)
            index=1;
        end
    end
    if(index==1)
        [~, Covariate,~]=pca(Covariate);
    end
    for i=1:1:size(Covariate,2)
        mu=mean(Covariate(:,i));
        sd=std(Covariate(:,i));
        Covariate(:,i)=(Covariate(:,i)-mu)/sd;
    end
    
    Change_point=NaN;
    
    [int_var,Model_type,int_MW_par]=initial_var(AMS,Covariate,best_fit,20);
    int_var=Variable_parameters1(int_var,Model_type,Change_point,AMS,int_MW_par,Covariate,best_fit);
    %% DEMC univariate Marginals evluation
    cha=5;
    evl=15000;
    bur=10000;
    sts=1;
    [para,Rhat,~,~]= demc_all(Model_type,int_var,AMS,evl,bur,sts,cha,Covariate,NaN,best_fit);
    x=1;
     while max(Rhat)>1.1||any(isnan(Rhat))%||accrate_PF<0.3
         x=x+1;
         if(x<=5)
            [para,Rhat,~,~]= demc_all(Model_type,int_var,AMS,evl,bur,sts,cha,Covariate,NaN,best_fit);
            disp(x);
         else
             break;
         end
     end
     [m,n]=size(Covariate);
     switch best_fit 
         case 'gev'
               
                % Taking exp of shape parameters;
                if(strcmp(Model_type(2,1),'Stationary'))
                    column=2;
                elseif(strcmp(Model_type(2,1),'Non-Stationary Trend'))
                    column=2:size(int_var,1);   
                end
                para1=para;
                % para(:,column)=exp(para(:,column)); %commented 
                PARA_FINAL(1:25005,1:size(para,2),com)=para;
                if(strcmp(Model_type(2,1),'Stationary') && strcmp(Model_type(3,1),'Stationary'))
                    shape=repmat(para(:,1),1,m);
                    scale_t=repmat(para(:,2),1,m);
                    mu_t=repmat(para(:,3),1,m);
                elseif(strcmp(Model_type(2,1),'Stationary') && strcmp(Model_type(3,1),'Non-Stationary Trend'))
                    shape=repmat(para(:,1),1,m);
                    scale_t = repmat(para(:,2),1,m);
                    mu_t = para(:,3:end-1)*Covariate' + para(:,end);
                elseif(strcmp(Model_type(2,1),'Non-Stationary Trend') && strcmp(Model_type(3,1),'Stationary'))
                    shape=repmat(para(:,1),1,m);
                    scale_t= para(:,2:n+1)*Covariate' + para(:,n+2);
                    mu_t = repmat(para(:,end),1,m);
                elseif(strcmp(Model_type(2,1),'Non-Stationary Trend') && strcmp(Model_type(3,1),'Non-Stationary Trend'))
                    shape=repmat(para(:,1),1,m);
                    scale_t= para(:,2:n+1)*Covariate' + para(:,n+2);
                    mu_t = para(:,n+3:end-1)*Covariate'+ para(:,end);
                end
                param(1,:)=median(shape);
                param(2,:)=median(scale_t);
                param(3,:)=median(mu_t);
                
                
                likelihood(1:size(AMS,1),1)=NaN;
                for s=1:1:size(AMS,1)
                    likelihood(s,1)=gevpdf(AMS(sts,1),param(1,s),exp(param(2,s)),param(3,s));
                end
         case 'gamma'
              % Taking exp of a and b parameters;
                para1=para;
                % para=exp(para); %commented
                PARA_FINAL(1:25005,1:size(para,2),com)=para;
                if(strcmp(Model_type(1,1),'Stationary') && strcmp(Model_type(2,1),'Stationary'))
                    shape=repmat(para(:,1),1,m);
                    scale_t=repmat(para(:,2),1,m);
                    % mu_t=repmat(para(:,3),1,m);
                elseif(strcmp(Model_type(1,1),'Stationary') && strcmp(Model_type(2,1),'Non-Stationary Trend'))
                    shape=repmat(para(:,1),1,m);
                    scale_t = para(:,2:end-1)*Covariate' + para(:,end);
                    % mu_t = para(:,4:end)*Covariate' + para(:,3);
                elseif(strcmp(Model_type(1,1),'Non-Stationary Trend') && strcmp(Model_type(2,1),'Stationary'))
                    shape= para(:,1:end-2)*Covariate' + para(:,end-1);
                    scale_t= repmat(para(:,end),1,m);
                    % mu_t = repmat(para(:,end),1,m);
                elseif(strcmp(Model_type(1,1),'Non-Stationary Trend') && strcmp(Model_type(2,1),'Non-Stationary Trend'))
                    shape=para(:,1:n)*Covariate' + para(:,n+1);
                    scale_t= para(:,n+2:end-1)*Covariate' + para(:,end);
                    % mu_t = para(:,n+4:end)*Covariate'+ para(:,n+3);
                end
                param(1,:)=median(shape);
                param(2,:)=median(scale_t);
                % param(3,:)=median(mu_t);
                
                
                likelihood(1:size(AMS,1),1)=NaN;
                for s=1:1:size(AMS,1)
                    likelihood(s,1)=gampdf(AMS(sts,1),exp(param(1,s)),exp(param(2,s)));
                end
         case 'gp'
              % Taking exp of shape parameters;
               if(strcmp(Model_type(2,1),'Stationary'))
                    column=2;
                elseif(strcmp(Model_type(2,1),'Non-Stationary Trend'))
                    column=2:size(int_var,1);   
                end
                para1=para;
                % para(:,column)=exp(para(:,column)); %commented
                PARA_FINAL(1:25005,1:size(para,2),com)=para;
                if(strcmp(Model_type(1,1),'Stationary') && strcmp(Model_type(2,1),'Stationary'))
                    shape=repmat(para(:,1),1,m);
                    scale_t=repmat(para(:,2),1,m);
                elseif(strcmp(Model_type(1,1),'Stationary') && strcmp(Model_type(2,1),'Non-Stationary Trend'))
                    shape=repmat(para(:,1),1,m);
                    scale_t = para(:,2:end-1)*Covariate' + para(:,end);
                elseif(strcmp(Model_type(1,1),'Non-Stationary Trend') && strcmp(Model_type(2,1),'Stationary'))
                    shape= para(:,1:end-2)*Covariate' + para(:,end-1);
                    scale_t= repmat(para(:,end),1,m);
                elseif(strcmp(Model_type(1,1),'Non-Stationary Trend') && strcmp(Model_type(2,1),'Non-Stationary Trend'))
                    shape=para(:,1:n)*Covariate' + para(:,n+1);
                    scale_t= para(:,n+2:end-1)*Covariate' + para(:,end);
                end
                param(1,:)=median(shape);
                param(2,:)=median(scale_t);
                likelihood(1:size(AMS,1),1)=NaN;
                for s=1:1:size(AMS,1)
                    likelihood(s,1)=gppdf(AMS(sts,1),param(1,s),exp(param(2,s)));
                end
         case 'lognormal'
              % Taking exp of SD parameters;
                if(strcmp(Model_type(2,1),'Stationary'))
                    column=2;
                elseif(strcmp(Model_type(2,1),'Non-Stationary Trend'))
                    column=2:size(int_var,1);   
                end
                 para1=para; %commented
                % para(:,column)=exp(para(:,column));
                PARA_FINAL(1:25005,1:size(para,2),com)=para;
                if(strcmp(Model_type(1,1),'Stationary') && strcmp(Model_type(2,1),'Stationary'))
                    shape=repmat(para(:,1),1,m);
                    scale_t=repmat(para(:,2),1,m);
                elseif(strcmp(Model_type(1,1),'Stationary') && strcmp(Model_type(2,1),'Non-Stationary Trend'))
                    shape=repmat(para(:,1),1,m);
                    scale_t = para(:,2:end-1)*Covariate' + para(:,end);
                elseif(strcmp(Model_type(1,1),'Non-Stationary Trend') && strcmp(Model_type(2,1),'Stationary'))
                    shape= para(:,1:end-2)*Covariate' + para(:,end-1);
                    scale_t= repmat(para(:,end),1,m);
                elseif(strcmp(Model_type(1,1),'Non-Stationary Trend') && strcmp(Model_type(2,1),'Non-Stationary Trend'))
                    shape=para(:,1:n)*Covariate' + para(:,n+1);
                    scale_t= para(:,n+2:end-1)*Covariate' + para(:,end);
                end
                param(1,:)=median(shape);
                param(2,:)=median(scale_t);
                likelihood(1:size(AMS,1),1)=NaN;
                for s=1:1:size(AMS,1)
                    likelihood(s,1)=lognpdf(AMS(sts,1),param(1,s),exp(param(2,s)));
                end
    
     end
        
             k1 = numel(para(1,:));       % Number of parameters in model 1
             L1 = nansum(log(likelihood));  % Log-likelihood of model 1
        
             %%do something to store the name and AIC together in an array
             [DIC3]=DIC(para1(bur+1:end,:),Covariate,AMS,Model_type,best_fit);
             DIC1(com,1)=DIC3;
             DIC2{com,1}=sprintf('Covariate %d',currentCombination);
             AIC1(com,1) = ( 2*k1 - 2*L1);         % AIC for model 1
             AIC2{com,1} = sprintf('Covariate %d',currentCombination);
             displayed=sprintf('Model %d',com);
             % disp(displayed)
                com=com+1;
              toc;
            end
        end
    end


else

[int_var,Model_type,int_MW_par]=initial_var(AMS,Covariate1,best_fit,20);
int_var=Variable_parameters1(int_var,Model_type,NaN,AMS,int_MW_par,Covariate1,best_fit);
%% DEMC univariate Marginals evluation
cha=5;
evl=15000;
bur=10000;
sts=1;
[para,Rhat,~,~]= demc_all(Model_type,int_var,AMS,evl,bur,sts,cha,Covariate1,NaN,best_fit);
x=1;
 while max(Rhat)>1.1||any(isnan(Rhat))%||accrate_PF<0.3
     x=x+1;
     if(x<=5)
        [para,Rhat,~,~]= demc_all(Model_type,int_var,AMS,evl,bur,sts,cha,Covariate1,NaN,best_fit);
        disp(x);
     else
         break;
     end
 end
 [m,n]=size(Covariate1);
 switch best_fit 
     case 'gev'
           
            % Taking exp of shape parameters;
            if(strcmp(Model_type(2,1),'Stationary'))
                column=2;
            elseif(strcmp(Model_type(2,1),'Non-Stationary Trend'))
                column=2:size(int_var,1);   
            end
            para1=para;
            % para(:,column)=exp(para(:,column)); %commented
            PARA_FINAL(1:25005,1:size(para,2),com)=para; 
     case 'gamma'
            para1=para;
            % para=exp(para); %commented
            PARA_FINAL(1:25005,1:size(para,2),com)=para;
            
     case 'gp'
           if(strcmp(Model_type(2,1),'Stationary'))
                column=2;
            elseif(strcmp(Model_type(2,1),'Non-Stationary Trend'))
                column=2:size(int_var,1);   
            end
            para1=para;
            % para(:,column)=exp(para(:,column)); %commented
            PARA_FINAL(1:25005,1:size(para,2),com)=para;
            
     case 'lognormal'
            if(strcmp(Model_type(2,1),'Stationary'))
                column=2;
            elseif(strcmp(Model_type(2,1),'Non-Stationary Trend'))
                column=2:size(int_var,1);   
            end
            para1=para;
            % para(:,column)=exp(para(:,column)); %commented
            PARA_FINAL(1:25005,1:size(para,2),com)=para;
            
 end

DIC1=DIC(para1(bur+1:end,:),Covariate1,AMS,Model_type,best_fit);
DIC2=sprintf('NaN');
AIC1 =NaN;       % AIC for model 1
AIC2 = sprintf('NaN')
end
warning('on', 'all');
end
