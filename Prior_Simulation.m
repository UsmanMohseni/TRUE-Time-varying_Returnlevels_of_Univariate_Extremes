%% Prior estimator for simulation part
function [p_x,x_k]=Prior_Simulation(int_var,~,~,Covariate,obs_int,Variance,~,Z_Z,Model_type,pr_dist_name,hy_p1,hy_p2)

if(strcmp(Model_type,'Stationary'))
    x_k_temp= Z_Z(:,1);
    x_k=repmat(x_k_temp,1,length(obs_int));
    % p_x_k_temp= normpdf(x_k_temp,0,Variance);
    p_x_k_temp=eval([pr_dist_name,'pdf','(','x_k_temp',',','hy_p1',',','hy_p2',')']);

    p_x=log(p_x_k_temp);
   
elseif(strcmp(Model_type,'Non-Stationary Trend'))
    k_s= Z_Z(:,1:size(int_var,1)-1);
    k_i= Z_Z(:,end);
    % p_ks=normpdf(k_s,0,Variance);
    % p_ki=normpdf(k_i,0,Variance);
    p_ks=eval([pr_dist_name,'pdf','(','k_s',',','hy_p1',',','hy_p2',')']);
    p_ki=eval([pr_dist_name,'pdf','(','k_i',',','hy_p1',',','hy_p2',')']);

    for it=1:length(obs_int)
        x_k(:,it)= sum(k_s*Covariate(it,:)',2)+k_i; %changed
    end
    p_x=sum(log(p_ks)+log(p_ki),2);
end
end
