function [para_PF,Rhat_PF,accrate_PF,mix_PF]=demc_all_Updated(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN,best_fit)
switch best_fit
    case 'gev'
[para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_gev(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN) ;
% [para_PD,Rhat_PD,accrate_PD,mix_PD]= demc_gev(Model_type2,int2_var,PD,evl,bur,sts,cha,Covariate(:,Best_Covariate_PD),NaN) ;
% [para_PV,Rhat_PV,accrate_PV,mix_PV]= demc_gev(Model_type3,int3_var,PV,evl,bur,sts,cha,Covariate(:,Best_Covariate_PV),NaN) ;

    case 'gp'
[para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_gp(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN) ;
% [para_PD,Rhat_PD,accrate_PD,mix_PD]= demc_gp(Model_type2,int2_var,PD,evl,bur,sts,cha,Covariate(:,Best_Covariate_PD),NaN) ;
% [para_PV,Rhat_PV,accrate_PV,mix_PV]= demc_gp(Model_type3,int3_var,PV,evl,bur,sts,cha,Covariate(:,Best_Covariate_PV),NaN) ;

    case 'gamma'
[para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_gam(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN) ;
% [para_PD,Rhat_PD,accrate_PD,mix_PD]= demc_gam(Model_type2,int2_var,PD,evl,bur,sts,cha,Covariate(:,Best_Covariate_PD),NaN) ;
% [para_PV,Rhat_PV,accrate_PV,mix_PV]= demc_gam(Model_type3,int3_var,PV,evl,bur,sts,cha,Covariate(:,Best_Covariate_PV),NaN) ;

    case 'lognormal'
[para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_lognormal(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN) ;
% [para_PD,Rhat_PD,accrate_PD,mix_PD]= demc_lognormal(Model_type2,int2_var,PD,evl,bur,sts,cha,Covariate(:,Best_Covariate_PD),NaN) ;
% [para_PV,Rhat_PV,accrate_PV,mix_PV]= demc_lognormal(Model_type3,int3_var,PV,evl,bur,sts,cha,Covariate(:,Best_Covariate_PV),NaN) ;

    case 'normal'
[para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_normal(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN) ;

    case 'weibull'
[para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_wbl(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN) ;

    case 'ev'
[para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_ev(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN) ;

    case 'logistic'
[para_PF,Rhat_PF,accrate_PF,mix_PF]= demc_logistic(Model_type1,int1_var,PF,evl,bur,sts,cha,Covariate,NaN) ;

end
end
