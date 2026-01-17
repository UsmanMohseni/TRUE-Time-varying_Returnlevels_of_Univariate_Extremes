%% Function for Prior estimation used in DE-MC
function [p_x, mix_x,x_k]=Prior_estimator(int_var,~,~,Covariate,obs_int,~,cha,Model_type,pr_dist_name,hy_p1,hy_p2)
if(strcmp(Model_type,'Stationary'))
    %  count=0;
    % a=abs(int_var(1,1));
    % if(a<1)
    % while(a<1)
    %     a=a*10;
    %     count=count+1;
    % end
    % b=10^-(count+1);
    % elseif a>1
    % while(a>1)
    % a=a/10;
    % count=count+1;
    % end
    % 
    % b=10^(count-3);
    % end
    b=0.01*int_var(1,1);  %% Here insted of 0.001 we have took the noise scale to be 1 percent of the value. It will resolve the problem of nan
    x_k_temp= int_var(1,1)+b*rand(cha,1);
    x_k=repmat(x_k_temp,1,length(obs_int));

   p_x_k_temp=eval([pr_dist_name,'pdf','(','x_k_temp',',','hy_p1',',','hy_p2',')']);
    % p_x_k_temp= normpdf(x_k_temp,0,Variance);
    p_x=log(p_x_k_temp);
    mix_x(1:cha,1,1)= x_k_temp;
 elseif (strcmp(Model_type,'Non-Stationary Trend'))
    %   count=0;
    % a=abs(int_var(1,1));
    % while(a<1)
    %     a=a*10;
    %     count=count+1;
    % end
    % b=10^-(count+1);

    for i=1:size(int_var,1)-1
         b=0.01*int_var(i,1);%% Here insted of 0.001 we have took the noise scale to be 1 percent of the value. It will resolve the problem of nan
         k_s(:,i)= repmat(int_var(i,1),cha,1)+b*rand(cha,1);
    end
    b=0.01*int_var(end,1);%% Here insted of 0.001 we have took the noise scale to be 1 percent of the value. It will resolve the problem of nan
    k_i= repmat(int_var(end,1),cha,1)+b*rand(cha,1);
       p_ks=eval([pr_dist_name,'pdf','(','k_s',',','hy_p1',',','hy_p2',')']);
       p_ki=eval([pr_dist_name,'pdf','(','k_i',',','hy_p1',',','hy_p2',')']);

    % p_ks=normpdf(k_s,0,Variance);
    % p_ki=normpdf(k_i,0,Variance);
    for it=1:length(obs_int)
        x_k(:,it)= sum(k_s*Covariate(it,:)',2)+k_i;%changed
    end
    p_x=sum(log(p_ks)+log(p_ki),2);
    mix_x(1:cha,:,1)= [k_s k_i]; %changed
 end

end