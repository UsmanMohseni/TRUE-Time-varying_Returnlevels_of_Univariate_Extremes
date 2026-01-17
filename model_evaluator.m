function [DIC1,DIC2]= model_evaluator(numVariables,allCombinations,Covariate1,AMS)

for variableIndex = 1:numVariables
    combinationsForVariable = allCombinations{variableIndex};
    
    % Loop through each combination
    for combinationIndex = 1:size(combinationsForVariable, 1)
        currentCombination = combinationsForVariable(combinationIndex, :);
        if(currentCombination(:,end)~=7)
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




%Standardised AMS
AMS1=AMS;
mu1=mean(AMS1(:,1));
sd1=std(AMS1(:,1));
for i=1:length(AMS1)
AMS_standard(i,1)=(AMS1(i)-mu1)/sd1;
end
AMS=AMS_standard;


%% Model Type: Non-Stationary with Time as a covariate
%%Model with only Time as Covariate
% Covariate=DATA(:,end);
 [p,q]=size(Covariate);

int_var(1:q+1,1:6)=NaN;
MW=20;
Change_point=NaN;
[int_MW_par]=Par_series_TSW('gev',AMS,MW);
int_MW_par(2,:)=log(int_MW_par(2,:)); % log is taken to ensure the positive values of scale parameter in GEV
for i=1:1:3
    [Z(i,1) Z(i,2) Z(i,3)]=M_K_test(int_MW_par(i,:));
   % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
end
Z(1,end)=0;
%col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)

int_var(1,[1,3,5])=0;
k=1;
for i=1:1:3
   if Z(i,3)~=0
       int_var(2:end,k)=0;  
   end
        
   k=k+2;
end

% Defining Model type
% Stationary
% Non-Stationary Trend
% Step Trend
% Step No-Trend
%% Model Type: STATIONARY
%%Stationry Model
% Covariate=ones(length(AMS),1);
% [p,q]=size(Covariate);
%  int_var(2,5)=NaN;
% int_var(2,3)=NaN;
Model_type={'-','-';'-','-';'-','-'};



if(isnan(int_var(1,2))&&isnan(int_var(2,2))&&isnan(int_var(2,1)))
    Model_type{1,1}='Stationary';
end

if(isnan(int_var(1,4))&&isnan(int_var(2,4))&&isnan(int_var(2,3)))
    Model_type{2,1}='Stationary';
elseif(~isnan(int_var(1,4))&&isnan(int_var(2,4))&&isnan(int_var(2,3)))
    Model_type{2,1}='Step No-Trend';
    Model_type{2,2}='Step No-Trend';
elseif(isnan(int_var(1,4))&&isnan(int_var(2,4))&&~isnan(int_var(2,3)))
    Model_type{2,1}='Non-Stationary Trend';
elseif(~isnan(int_var(1,4))&&~isnan(int_var(2,4))&&isnan(int_var(2,3)))
    Model_type{2,1}='Step No-Trend';
    Model_type{2,2}='Step Trend';
elseif(~isnan(int_var(1,4))&&isnan(int_var(2,4))&&~isnan(int_var(2,3)))
    Model_type{2,1}='Step Trend';
    Model_type{2,2}='Step No-Trend';
else
    Model_type{2,1}='Step Trend';
    Model_type{2,2}='Step Trend';
end

if(isnan(int_var(1,6))&&isnan(int_var(2,6))&&isnan(int_var(2,5)))
    Model_type{3,1}='Stationary';
elseif(~isnan(int_var(1,6))&&isnan(int_var(2,6))&&isnan(int_var(2,5)))
    Model_type{3,1}='Step No-Trend';
    Model_type{3,2}='Step No-Trend';
elseif(isnan(int_var(1,6))&&isnan(int_var(2,6))&&~isnan(int_var(2,5)))
    Model_type{3,1}='Non-Stationary Trend';
elseif(~isnan(int_var(1,6))&&~isnan(int_var(2,6))&&isnan(int_var(2,5)))
    Model_type{3,1}='Step No-Trend';
    Model_type{3,2}='Step Trend';
elseif(~isnan(int_var(1,6))&&isnan(int_var(2,6))&&~isnan(int_var(2,5)))
    Model_type{3,1}='Step Trend';
    Model_type{3,2}='Step No-Trend';
else
    Model_type{3,1}='Step Trend';
    Model_type{3,2}='Step Trend';
end

%Preparing the variable parameters table
Sta_para=gevfit(AMS);
int_var(1,1)=Sta_para(1,1);
int_var(2,1)=NaN;
int_var(:,2)=NaN;




[int_var]=Variable_parameters(int_var,Model_type,Change_point,AMS,int_MW_par,Covariate);

%DE-MC Iterations
evl=7000;
cha=5;
bur=2000;
sts=1;

[para,Rhat,accrate,mix]= demc_gev(Model_type,int_var,AMS,evl,bur,sts,cha,Covariate,Change_point) ;
count=0;
while max(Rhat)>1.1||any(isnan(Rhat))
   count=count+1
   if(count<=5)
    [para,Rhat,accrate,mix]= demc_gev(Model_type,int_var,AMS,evl,bur,sts,cha,Covariate,Change_point) ;
   else

       break;
   end
end

[m,n]=size(Covariate);
% Taking exp of shape parameters;
if(strcmp(Model_type(2,1),'Stationary'))
    column=[2];
elseif(strcmp(Model_type(2,1),'Non-Stationary Trend'))
    column=[2:size(int_var,1)];   
end
para1=para;
 para(:,column)=exp(para(:,column));

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
%% changes
param(1,:)=median(shape);
param(2,:)=median(scale_t);
param(3,:)=median(mu_t);

%%


k1 = numel(para(1,:));       % Number of parameters in model 1
for s=1:1:size(AMS)
likelihood(s,1)=gevpdf(AMS(sts,1),param(1,s),param(2,s),param(3,s));
end
L1 = nansum(log(likelihood));  % Log-likelihood of model 1

% %do something to store the name and AIC together in an array
[DIC3]=DIC(para1(10001:end,:),Covariate,AMS,Model_type);
DIC1(com,1)=DIC3;
DIC2{com,1}=sprintf('Covariate %d',currentCombination);
AIC1(com,1) = ( 2*k1 - 2*L1);         % AIC for model 1
AIC2{com,1} = sprintf('Covariate %d',currentCombination);
displayed=sprintf('Model %d',com);
disp(displayed)
com=com+1;
       
    end
    end
end
end
