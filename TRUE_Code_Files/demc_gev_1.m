%%  Function for DE-MC
function [para,Rhat,accrate,mix]= demc_gev_1(Model_type,int_var,AMS,evl,bur,sts,cha,Covariate,int_cp_index) 
% rng(rng_index)
% function [para,Rhat,accrate,mix]= demc_gev(Model_type,int_var,AMS,evl,bur,sts,cha,Covariate,int_cp_index,rng_index) 
% rng(rng_index)
%finding the Number of parameters (dim)
[m,n]=size(int_var);
dim=0;
for i=1:1:m
    for j=1:1:n
        if(~isnan(int_var(i,j)))
            dim=dim+1;
        end
    end
end


% DE-MC simulation part
accept= zeros(evl,1);
reject= zeros(evl,1);
mix= zeros(cha,dim+1,evl);
z_z= zeros(cha,dim);
obs_int= repmat(AMS',cha,1);

% % % Preparing Initial Chains
% % [p_x_k,mix_k,x_k]=Prior_estimator(int_var(:,1:2),int_cp_index,AMS,Covariate,obs_int,0.3,cha,Model_type{1,1},rng_index);
% % [ p_x_sigma, mix_x_sig,x_sigma]=Prior_estimator(int_var(:,3:4),int_cp_index,AMS,Covariate,obs_int,100,cha,Model_type{2,1},rng_index);
% % [ p_x_mu, mix_x_mu,x_mu]=Prior_estimator(int_var(:,5:6),int_cp_index,AMS,Covariate,obs_int,100,cha,Model_type{3,1},rng_index);
% Preparing Initial Chains
[p_x_k,mix_k,x_k]=Prior_estimator(int_var(:,1:2),int_cp_index,AMS,Covariate,obs_int,0.3,cha,Model_type{1,1});
[ p_x_sigma, mix_x_sig,x_sigma]=Prior_estimator(int_var(:,3:4),int_cp_index,AMS,Covariate,obs_int,100,cha,Model_type{2,1});
[ p_x_mu, mix_x_mu,x_mu]=Prior_estimator(int_var(:,5:6),int_cp_index,AMS,Covariate,obs_int,100,cha,Model_type{3,1});

non_xx= (log(gevpdf(obs_int,x_k,exp(x_sigma),x_mu)));
likeli_joint= sum(non_xx,2);
p_post= likeli_joint+p_x_sigma+p_x_k+p_x_mu;
mix(1:cha,1:dim,1)= [mix_k mix_x_sig mix_x_mu];
mix(1:cha,end,1)= p_post;
JR(:,1)= 2.38./sqrt(2*(1:dim)');

% Monte-Carlo Simulation
for i= 1:evl
    
    [~,ttt]= sort(rand(cha-1,cha));
    D = rand(cha,dim);
    for jj= 1:cha
        ii= ones(cha,1); ii(jj)= 0; 
        idx= find(ii> 0);
        rr= idx(ttt(1:2,jj));
        [c]= find(D(jj,1:dim)> 0);
        Dim= size(c,2);
        if isempty(c), c= randperm(dim); c= c(1);
        end
        if rand< 0.8
            g= JR(Dim,1);
        else
            g= 1;
        end
        delta(jj,c)= g*((mix(rr(1),1:dim,i)-mix(rr(2),1:dim,i)));
        if (sum(delta(jj,1:dim).^2,2)==0)
            [RR,P]= chol(cov(mix(:,1:dim,i))+1e-5*eye(dim));
            if P==0
                R= (2.38/sqrt(dim))*RR;
                delta(jj,1:dim)= randn(1,dim)*R;
            end
        end
    end
    z_z(:,1:dim)= mix(:,1:dim,i)+delta(:,1:dim);
    p_post= mix(:,end,i);

    if(strcmp(Model_type(2,1),'Stationary')&& strcmp(Model_type(3,1),'Stationary'))

    [p_z_k,z_k]=Prior_Simulation(int_var(:,1:2),int_cp_index,AMS,Covariate,obs_int,0.3,cha,z_z(:,1),Model_type{1,1});
    [p_z_sigma,z_sigma]=Prior_Simulation(int_var(:,3:4),int_cp_index,AMS,Covariate,obs_int,100,cha,z_z(:,2),Model_type{2,1});
    [p_z_mu,z_mu]=Prior_Simulation(int_var(:,5:6),int_cp_index,AMS,Covariate,obs_int,100,cha,z_z(:,3),Model_type{3,1});
    elseif (strcmp(Model_type(2,1),'Stationary')&& strcmp(Model_type(3,1),'Non-Stationary Trend'))
    [p_z_k,z_k]=Prior_Simulation(int_var(:,1:2),int_cp_index,AMS,Covariate,obs_int,0.3,cha,z_z(:,1),Model_type{1,1});
    [p_z_sigma,z_sigma]=Prior_Simulation(int_var(:,3:4),int_cp_index,AMS,Covariate,obs_int,100,cha,z_z(:,2),Model_type{2,1});
    [p_z_mu,z_mu]=Prior_Simulation(int_var(:,5:6),int_cp_index,AMS,Covariate,obs_int,100,cha,z_z(:,3:end),Model_type{3,1});
    elseif (strcmp(Model_type(2,1),'Non-Stationary Trend')&& strcmp(Model_type(3,1),'Stationary'))
    [p_z_k,z_k]=Prior_Simulation(int_var(:,1:2),int_cp_index,AMS,Covariate,obs_int,0.3,cha,z_z(:,1),Model_type{1,1});
    [p_z_sigma,z_sigma]=Prior_Simulation(int_var(:,3:4),int_cp_index,AMS,Covariate,obs_int,100,cha,z_z(:,2:m+1),Model_type{2,1});
    [p_z_mu,z_mu]=Prior_Simulation(int_var(:,5:6),int_cp_index,AMS,Covariate,obs_int,100,cha,z_z(:,m+2:end),Model_type{3,1});
    elseif(strcmp(Model_type(2,1),'Non-Stationary Trend')&& strcmp(Model_type(3,1),'Non-Stationary Trend'))
    [p_z_k,z_k]=Prior_Simulation(int_var(:,1:2),int_cp_index,AMS,Covariate,obs_int,0.3,cha,z_z(:,1),Model_type{1,1});
    [p_z_sigma,z_sigma]=Prior_Simulation(int_var(:,3:4),int_cp_index,AMS,Covariate,obs_int,100,cha,z_z(:,2:m+1),Model_type{2,1});
    [p_z_mu,z_mu]=Prior_Simulation(int_var(:,5:6),int_cp_index,AMS,Covariate,obs_int,100,cha,z_z(:,m+2:end),Model_type{3,1});
    end
    non_zxx= (log(gevpdf(obs_int,z_k,exp(z_sigma),z_mu)));
    p_likeli_joint= sum(non_zxx,2);
    pz_post=p_likeli_joint+p_z_sigma+p_z_k+p_z_mu;
       ratio= pz_post-p_post;
    alfa= (min(exp(ratio),1));
    r= rand(cha,1);
    cidx= find(alfa>r);
    mix(cidx,1:dim,i+1)= z_z(cidx,1:dim);
    mix(cidx,end,i+1)= pz_post(cidx);
    ridx= find(alfa<=r);
    mix(ridx,1:dim,i+1)= mix(ridx,1:dim,i);
    mix(ridx,end,i+1)= mix(ridx,end,i);
    accept(i,1)= size(cidx,1);
    reject(i,1)= size(ridx,1);
     % disp(i)
end


for i=1:1:dim
    P=mix(:,i,bur:evl);
    P1= reshape(P,1,cha*(evl-bur+1),1);
    para(:,i)=P1(:,sts:end)';
end
sur='nonsta';
Rhat= Rconver(1,evl,bur,cha,mix,sts,sur,dim);
accrate= sum(accept)./sum((accept+reject));
%disp(accrate)

end

