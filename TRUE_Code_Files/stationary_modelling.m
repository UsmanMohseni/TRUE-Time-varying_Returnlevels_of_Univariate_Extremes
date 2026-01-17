function [param_PF_shape, param_PF_scale_t, param_PF_mu_t,...
    param_PD_shape, param_PD_scale_t, param_PD_mu_t,...
    param_PV_shape, param_PV_scale_t, param_PV_mu_t,...
    teta_PF_PD,teta_PF_PV,teta_PV_PD]=stationary_modelling(DATA)
    %% Standardizing the DATA
[o w]=size(DATA);
DATA_Standardized(1:o,1:w)=NaN;
DATA_Standardized(:,[1,end])=DATA(:,[1,end]);

mean_PF=mean(DATA(:,2));
sd_PF=sqrt(var(DATA(:,2)));

mean_PD=mean(DATA(:,3));
sd_PD=sqrt(var(DATA(:,3)));

mean_PV=mean(DATA(:,4));
sd_PV=sqrt(var(DATA(:,4)));

DATA_Standardized(:,2)=(DATA(:,2)-mean_PF)./sd_PF;
DATA_Standardized(:,3)=(DATA(:,3)-mean_PD)./sd_PD;
DATA_Standardized(:,4)=(DATA(:,4)-mean_PV)./sd_PV;

for i=5:w-1
    avg=mean(DATA(:,i));
    SD=sqrt(var(DATA(:,i)));
    DATA_Standardized(:,i)=(DATA(:,i)-avg)./SD;
end

%% Seperating all data
Years=DATA_Standardized(:,1);
PF=DATA_Standardized(:,2);
PD=DATA_Standardized(:,3);
PV=DATA_Standardized(:,4);
Covariate=DATA_Standardized(:,5:end-1);
Range=DATA_Standardized(:,1);

%% Variables parameter series
Covariate=ones(length(PF),1);
[p,q]=size(Covariate);
int1_var(1:q+1,1:6)=NaN;
MW=20;
Change_point=NaN;
[int1_MW_par]=Par_series_TSW('gev',PF,MW);
int1_MW_par(2,:)=log(int1_MW_par(2,:)); % log is taken to ensure the positive values of scale parameter in GEV
for i=1:1:3
    [Z1(i,1) Z1(i,2) Z1(i,3)]=M_K_test(int1_MW_par(i,:));
   % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
end
Z1(1,end)=0;
%col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)

int1_var(1,[1,3,5])=0;
k=1;
for i=1:1:3
   if Z1(i,3)~=0
       int1_var(2:end,k)=0;  
   end

   k=k+2;
end

% Defining Model type
% Stationary
% Non-Stationary Trend
% Step Trend
% Step No-Trend
%% Model Type: STATIONARY
%Stationry Model
int1_var(2,5)=NaN;
int1_var(2,3)=NaN;
Model_type1={'-','-';'-','-';'-','-'};


if(isnan(int1_var(1,2))&&isnan(int1_var(2,2))&&isnan(int1_var(2,1)))
    Model_type1{1,1}='Stationary';
end

if(isnan(int1_var(1,4))&&isnan(int1_var(2,4))&&isnan(int1_var(2,3)))
    Model_type1{2,1}='Stationary';
elseif(~isnan(int1_var(1,4))&&isnan(int1_var(2,4))&&isnan(int1_var(2,3)))
    Model_type1{2,1}='Step No-Trend';
    Model_type1{2,2}='Step No-Trend';
elseif(isnan(int1_var(1,4))&&isnan(int1_var(2,4))&&~isnan(int1_var(2,3)))
    Model_type1{2,1}='Non-Stationary Trend';
elseif(~isnan(int1_var(1,4))&&~isnan(int1_var(2,4))&&isnan(int1_var(2,3)))
    Model_type1{2,1}='Step No-Trend';
    Model_type1{2,2}='Step Trend';
elseif(~isnan(int1_var(1,4))&&isnan(int1_var(2,4))&&~isnan(int1_var(2,3)))
    Model_type1{2,1}='Step Trend';
    Model_type1{2,2}='Step No-Trend';
else
    Model_type1{2,1}='Step Trend';
    Model_type1{2,2}='Step Trend';
end

if(isnan(int1_var(1,6))&&isnan(int1_var(2,6))&&isnan(int1_var(2,5)))
    Model_type1{3,1}='Stationary';
elseif(~isnan(int1_var(1,6))&&isnan(int1_var(2,6))&&isnan(int1_var(2,5)))
    Model_type1{3,1}='Step No-Trend';
    Model_type1{3,2}='Step No-Trend';
elseif(isnan(int1_var(1,6))&&isnan(int1_var(2,6))&&~isnan(int1_var(2,5)))
    Model_type1{3,1}='Non-Stationary Trend';
elseif(~isnan(int1_var(1,6))&&~isnan(int1_var(2,6))&&isnan(int1_var(2,5)))
    Model_type1{3,1}='Step No-Trend';
    Model_type1{3,2}='Step Trend';
elseif(~isnan(int1_var(1,6))&&isnan(int1_var(2,6))&&~isnan(int1_var(2,5)))
    Model_type1{3,1}='Step Trend';
    Model_type1{3,2}='Step No-Trend';
else
    Model_type1{3,1}='Step Trend';
    Model_type1{3,2}='Step Trend';
end

%Preparing the variable parameters table
Sta_para=gevfit(PF);
int1_var(1,1)=Sta_para(1,1);
int1_var(2,1)=NaN;
int1_var(:,2)=NaN;


int2_var(1:q+1,1:6)=NaN;
MW=20;
Change_point=NaN;
[int2_MW_par]=Par_series_TSW('gev',PD,MW);
int2_MW_par(2,:)=log(int2_MW_par(2,:)); % log is taken to ensure the positive values of scale parameter in GEV
for i=1:1:3
    [Z2(i,1) Z2(i,2) Z2(i,3)]=M_K_test(int2_MW_par(i,:));
   % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
end
Z2(1,end)=0;
%col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)

int2_var(1,[1,3,5])=0;
k=1;
for i=1:1:3
   if Z2(i,3)~=0
       int2_var(2:end,k)=0;  
   end

   k=k+2;
end

% Defining Model type
% Stationary
% Non-Stationary Trend
% Step Trend
% Step No-Trend
%% Model Type: STATIONARY
%Stationry Model

int2_var(2,5)=NaN;
int2_var(2,3)=NaN;
Model_type2={'-','-';'-','-';'-','-'};


if(isnan(int2_var(1,2))&&isnan(int2_var(2,2))&&isnan(int2_var(2,1)))
    Model_type2{1,1}='Stationary';
end

if(isnan(int2_var(1,4))&&isnan(int2_var(2,4))&&isnan(int2_var(2,3)))
    Model_type2{2,1}='Stationary';
elseif(~isnan(int2_var(1,4))&&isnan(int2_var(2,4))&&isnan(int2_var(2,3)))
    Model_type2{2,1}='Step No-Trend';
    Model_type2{2,2}='Step No-Trend';
elseif(isnan(int2_var(1,4))&&isnan(int2_var(2,4))&&~isnan(int2_var(2,3)))
    Model_type2{2,1}='Non-Stationary Trend';
elseif(~isnan(int2_var(1,4))&&~isnan(int2_var(2,4))&&isnan(int2_var(2,3)))
    Model_type2{2,1}='Step No-Trend';
    Model_type2{2,2}='Step Trend';
elseif(~isnan(int2_var(1,4))&&isnan(int2_var(2,4))&&~isnan(int2_var(2,3)))
    Model_type2{2,1}='Step Trend';
    Model_type2{2,2}='Step No-Trend';
else
    Model_type2{2,1}='Step Trend';
    Model_type2{2,2}='Step Trend';
end

if(isnan(int2_var(1,6))&&isnan(int2_var(2,6))&&isnan(int2_var(2,5)))
    Model_type2{3,1}='Stationary';
elseif(~isnan(int2_var(1,6))&&isnan(int2_var(2,6))&&isnan(int2_var(2,5)))
    Model_type2{3,1}='Step No-Trend';
    Model_type2{3,2}='Step No-Trend';
elseif(isnan(int2_var(1,6))&&isnan(int2_var(2,6))&&~isnan(int2_var(2,5)))
    Model_type2{3,1}='Non-Stationary Trend';
elseif(~isnan(int2_var(1,6))&&~isnan(int2_var(2,6))&&isnan(int2_var(2,5)))
    Model_type2{3,1}='Step No-Trend';
    Model_type2{3,2}='Step Trend';
elseif(~isnan(int2_var(1,6))&&isnan(int2_var(2,6))&&~isnan(int2_var(2,5)))
    Model_type2{3,1}='Step Trend';
    Model_type2{3,2}='Step No-Trend';
else
    Model_type2{3,1}='Step Trend';
    Model_type2{3,2}='Step Trend';
end

%Preparing the variable parameters table
Sta_para=gevfit(PD);
int2_var(1,1)=Sta_para(1,1);
int2_var(2,1)=NaN;
int2_var(:,2)=NaN;

int3_var(1:q+1,1:6)=NaN;
MW=20;
Change_point=NaN;
[int3_MW_par]=Par_series_TSW('gev',PV,MW);
int3_MW_par(2,:)=log(int3_MW_par(2,:)); % log is taken to ensure the positive values of scale parameter in GEV
for i=1:1:3
    [Z3(i,1) Z3(i,2) Z3(i,3)]=M_K_test(int3_MW_par(i,:));
   % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
end
Z3(1,end)=0;
%col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)

int3_var(1,[1,3,5])=0;
k=1;
for i=1:1:3
   if Z3(i,3)~=0
       int3_var(2:end,k)=0;  
   end

   k=k+2;
end

% Defining Model type
% Stationary
% Non-Stationary Trend
% Step Trend
% Step No-Trend
%% Model Type: STATIONARY
%Stationry Model

int3_var(2,5)=NaN;
int3_var(2,3)=NaN;
Model_type3={'-','-';'-','-';'-','-'};


if(isnan(int3_var(1,2))&&isnan(int3_var(2,2))&&isnan(int3_var(2,1)))
    Model_type3{1,1}='Stationary';
end

if(isnan(int3_var(1,4))&&isnan(int3_var(2,4))&&isnan(int3_var(2,3)))
    Model_type3{2,1}='Stationary';
elseif(~isnan(int3_var(1,4))&&isnan(int3_var(2,4))&&isnan(int3_var(2,3)))
    Model_type3{2,1}='Step No-Trend';
    Model_type3{2,2}='Step No-Trend';
elseif(isnan(int3_var(1,4))&&isnan(int3_var(2,4))&&~isnan(int3_var(2,3)))
    Model_type3{2,1}='Non-Stationary Trend';
elseif(~isnan(int3_var(1,4))&&~isnan(int3_var(2,4))&&isnan(int3_var(2,3)))
    Model_type3{2,1}='Step No-Trend';
    Model_type3{2,2}='Step Trend';
elseif(~isnan(int3_var(1,4))&&isnan(int3_var(2,4))&&~isnan(int3_var(2,3)))
    Model_type3{2,1}='Step Trend';
    Model_type3{2,2}='Step No-Trend';
else
    Model_type3{2,1}='Step Trend';
    Model_type3{2,2}='Step Trend';
end

if(isnan(int3_var(1,6))&&isnan(int3_var(2,6))&&isnan(int3_var(2,5)))
    Model_type3{3,1}='Stationary';
elseif(~isnan(int3_var(1,6))&&isnan(int3_var(2,6))&&isnan(int3_var(2,5)))
    Model_type3{3,1}='Step No-Trend';
    Model_type3{3,2}='Step No-Trend';
elseif(isnan(int3_var(1,6))&&isnan(int3_var(2,6))&&~isnan(int3_var(2,5)))
    Model_type3{3,1}='Non-Stationary Trend';
elseif(~isnan(int3_var(1,6))&&~isnan(int3_var(2,6))&&isnan(int3_var(2,5)))
    Model_type3{3,1}='Step No-Trend';
    Model_type3{3,2}='Step Trend';
elseif(~isnan(int3_var(1,6))&&isnan(int3_var(2,6))&&~isnan(int3_var(2,5)))
    Model_type3{3,1}='Step Trend';
    Model_type3{3,2}='Step No-Trend';
else
    Model_type3{3,1}='Step Trend';
    Model_type3{3,2}='Step Trend';
end

%Preparing the variable parameters table
Sta_para=gevfit(PV);
int3_var(1,1)=Sta_para(1,1);
int3_var(2,1)=NaN;
int3_var(:,2)=NaN;
%% Copula parameter series

PD_copula_MW_par=copula_Par_series_TSW('gev',[PF PD],20);
PV_copula_MW_par=copula_Par_series_TSW('gev',[PF PV],20);
DV_copula_MW_par=copula_Par_series_TSW('gev',[PV PD],20);

% PDV_copula_MW_par=copula_Par_series_TSW('gev',[PF PD PV],38);

% copula_MW_par description
% Row 1: Copula parameters for Peak Flow and Peak Duration
% Row 2: Copula parameters for Peak Flow and Peak Volumne
% Row 3: Copula parameters for Peak Duration and Peak Volumne
copula_MW_par=[PD_copula_MW_par;PV_copula_MW_par;DV_copula_MW_par];


% copula_int_var description
% col 1 & col 2: PF and PD copula
% col 3 & col 4: PF and PV copula
% col 5 & col 6: PV and PD copula
copula_int_var(1:q+1,1:6)=NaN;
for i=1:1:3
    [Z_copula(i,1), Z_copula(i,2), Z_copula(i,3)]=M_K_test(copula_MW_par(i,:));
   % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
end
%col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
copula_int_var(1,[1,3,5])=0;
k=1;
for i=1:1:3
   if Z_copula(i,3)~=0
       copula_int_var(2:end,k)=0;  
   end

   k=k+2;
end
% Defining Model type
% Stationary
% Non-Stationary Trend
% Step Trend
% Step No-Trend
%% Model Type: STATIONARY
%Stationry Model

copula_int_var(2,5)=NaN;
copula_int_var(2,3)=NaN;
copula_int_var(2,1)=NaN;

% copula_Model_type Description:
% Row 1: PF and PD copula
% Row 2: PF and PV copula
% Row 3: PV and PD copula
copula_Model_type={'-','-';'-','-';'-','-'};


if(isnan(copula_int_var(1,2))&&isnan(copula_int_var(2,2))&&isnan(copula_int_var(2,1)))
    copula_Model_type{1,1}='Stationary';
elseif(~isnan(copula_int_var(1,2))&&isnan(copula_int_var(2,2))&&isnan(copula_int_var(2,1)))
    copula_Model_type{1,1}='Step No-Trend';
    copula_Model_type{1,2}='Step No-Trend';
elseif(isnan(copula_int_var(1,2))&&isnan(copula_int_var(2,2))&&~isnan(copula_int_var(2,1)))
    copula_Model_type{1,1}='Non-Stationary Trend';
elseif(~isnan(copula_int_var(1,2))&&~isnan(copula_int_var(2,2))&&isnan(copula_int_var(2,1)))
    copula_Model_type{1,1}='Step No-Trend';
    copula_Model_type{1,2}='Step Trend';
elseif(~isnan(copula_int_var(1,2))&&isnan(copula_int_var(2,2))&&~isnan(copula_int_var(2,1)))
    copula_Model_type{1,1}='Step Trend';
    copula_Model_type{1,2}='Step No-Trend';
else
    copula_Model_type{1,1}='Step Trend';
    copula_Model_type{1,2}='Step Trend';
end

if(isnan(copula_int_var(1,4))&&isnan(copula_int_var(2,4))&&isnan(copula_int_var(2,3)))
    copula_Model_type{2,1}='Stationary';
elseif(~isnan(copula_int_var(1,4))&&isnan(copula_int_var(2,4))&&isnan(copula_int_var(2,3)))
    copula_Model_type{2,1}='Step No-Trend';
    copula_Model_type{2,2}='Step No-Trend';
elseif(isnan(copula_int_var(1,4))&&isnan(copula_int_var(2,4))&&~isnan(copula_int_var(2,3)))
    copula_Model_type{2,1}='Non-Stationary Trend';
elseif(~isnan(copula_int_var(1,4))&&~isnan(copula_int_var(2,4))&&isnan(copula_int_var(2,3)))
    copula_Model_type{2,1}='Step No-Trend';
    copula_Model_type{2,2}='Step Trend';
elseif(~isnan(copula_int_var(1,4))&&isnan(copula_int_var(2,4))&&~isnan(copula_int_var(2,3)))
    copula_Model_type{2,1}='Step Trend';
    copula_Model_type{2,2}='Step No-Trend';
else
    copula_Model_type{2,1}='Step Trend';
    copula_Model_type{2,2}='Step Trend';
end

if(isnan(copula_int_var(1,6))&&isnan(copula_int_var(2,6))&&isnan(copula_int_var(2,5)))
    copula_Model_type{3,1}='Stationary';
elseif(~isnan(copula_int_var(1,6))&&isnan(copula_int_var(2,6))&&isnan(copula_int_var(2,5)))
    copula_Model_type{3,1}='Step No-Trend';
    copula_Model_type{3,2}='Step No-Trend';
elseif(isnan(copula_int_var(1,6))&&isnan(copula_int_var(2,6))&&~isnan(copula_int_var(2,5)))
    copula_Model_type{3,1}='Non-Stationary Trend';
elseif(~isnan(copula_int_var(1,6))&&~isnan(copula_int_var(2,6))&&isnan(copula_int_var(2,5)))
    copula_Model_type{3,1}='Step No-Trend';
    copula_Model_type{3,2}='Step Trend';
elseif(~isnan(copula_int_var(1,6))&&isnan(copula_int_var(2,6))&&~isnan(copula_int_var(2,5)))
    copula_Model_type{3,1}='Step Trend';
    copula_Model_type{3,2}='Step No-Trend';
else
    copula_Model_type{3,1}='Step Trend';
    copula_Model_type{3,2}='Step Trend';
end

%% DEMC parameter evaluation
cha=5;
evl=10000;
bur=5000;
sts=1;
int1_var=Variable_parameters(int1_var,Model_type1,Change_point,PF,int1_MW_par,Covariate);
int2_var=Variable_parameters(int2_var,Model_type2,Change_point,PD,int2_MW_par,Covariate);
int3_var=Variable_parameters(int3_var,Model_type3,Change_point,PV,int3_MW_par,Covariate);


[para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_gev(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN) ;
count=0;
[para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_gev(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN) ;

while max(Rhat_PF)>1.1||any(isnan(Rhat_PF))||accrate_PF<0.3
    clc
    [para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_gev(Model_type1,int1_var,PD,evl,bur,sts,cha,Covariate,NaN) ;
end


[para_PD,Rhat_PD,accrate_PD,mix_PD]= demc_gev(Model_type2,int2_var,PD,evl,bur,sts,cha,Covariate,NaN) ;

while max(Rhat_PD)>1.1||any(isnan(Rhat_PD))||accrate_PD<0.3
    clc
[para_PD,Rhat_PD,accrate_PD,mix_PD]= demc_gev(Model_type2,int2_var,PD,evl,bur,sts,cha,Covariate,NaN) ;
end


[para_PV,Rhat_PV,accrate_PV,mix_PV]= demc_gev(Model_type3,int3_var,PV,evl,bur,sts,cha,Covariate,NaN) ;

while max(Rhat_PV)>1.1||any(isnan(Rhat_PV))||accrate_PV<0.3
    clc
[para_PV,Rhat_PV,accrate_PV,mix_PV]= demc_gev(Model_type3,int3_var,PV,evl,bur,sts,cha,Covariate,NaN) ;
end

%% copula part
[PF_k, PF_sigma, PF_mu]=coupula_prior_estimator(Covariate,PF,Model_type1,para_PF,5);
[PD_k, PD_sigma, PD_mu]=coupula_prior_estimator(Covariate,PD,Model_type2,para_PD,5);
[PV_k, PV_sigma, PV_mu]=coupula_prior_estimator(Covariate,PV,Model_type3,para_PV,5);


[copula_int_var(:,1:2)]=variable_parameters_copula(copula_int_var(:,1:2),copula_Model_type(1,:),Change_point,PF,PV,copula_MW_par(1,:),Covariate);
[copula_int_var(:,3:4)]=variable_parameters_copula(copula_int_var(:,3:4),copula_Model_type(2,:),Change_point,PF,PD,copula_MW_par(2,:),Covariate);
[copula_int_var(:,5:6)]=variable_parameters_copula(copula_int_var(:,5:6),copula_Model_type(3,:),Change_point,PD,PV,copula_MW_par(3,:),Covariate);


% rng("default")
[copula_para_PF_PD,copula_Rhat_PF_PD,copula_accrate_PF_PD,copula_mix_PF_PD]=copula_demc_gev1(copula_Model_type{1,1},copula_int_var(:,1:2),PF,PD,PF_k,PF_sigma,PF_mu,PD_k,PD_sigma,PD_mu,cha,evl,sts,bur,Covariate);
[copula_para_PF_PV,copula_Rhat_PF_PV,copula_accrate_PF_PV,copula_mix_PF_PV]=copula_demc_gev1(copula_Model_type{2,1},copula_int_var(:,3:4),PF,PV,PF_k,PF_sigma,PF_mu,PV_k,PV_sigma,PV_mu,cha,evl,sts,bur,Covariate);
[copula_para_PV_PD,copula_Rhat_PV_PD,copula_accrate_PV_PD,copula_mix_PV_PD]=copula_demc_gev1(copula_Model_type{3,1},copula_int_var(:,5:6),PV,PD,PV_k,PV_sigma,PV_mu,PD_k,PD_sigma,PD_mu,cha,evl,sts,bur,Covariate);

% Calculating Time varying copula parameters for all three combinations
teta_PF_PD=copula_parameters_NS(copula_Model_type(1,1),copula_para_PF_PD,Covariate);
teta_PF_PV=copula_parameters_NS(copula_Model_type(2,1),copula_para_PF_PV,Covariate);
teta_PV_PD=copula_parameters_NS(copula_Model_type(3,1),copula_para_PV_PD,Covariate);

% Calculating Time varying variable parameters for all flood characterstics
[param_PF_shape, param_PF_scale_t, param_PF_mu_t]=variable_parameters_NS(Model_type1,para_PF,Covariate);
[param_PD_shape, param_PD_scale_t, param_PD_mu_t]=variable_parameters_NS(Model_type2,para_PD,Covariate);
[param_PV_shape, param_PV_scale_t, param_PV_mu_t]=variable_parameters_NS(Model_type3,para_PV,Covariate);

    
end