function [DIC1]=DIC_Updated(para,Covariate,AMS,Model_type,best_fit)
%% calculating L Value
AMS1=AMS';
AMS=AMS1;
np=size(Covariate,2);
n=np;
avg_para=mean(para);


 switch best_fit 
     case 'gev'
            if strcmp(Model_type{2,1},'Stationary')&& strcmp(Model_type{3,1},'Stationary')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=repmat(avg_para(1,2),1,length(AMS));
                x_mu=repmat(avg_para(1,3),1,length(AMS));
            elseif strcmp(Model_type{2,1},'Stationary')&&strcmp(Model_type{3,1},'Non-Stationary Trend')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=repmat(avg_para(1,2),1,length(AMS));
                x_mu=avg_para(1,3:end-1)*Covariate'+avg_para(1,end);
            elseif strcmp(Model_type{2,1},'Non-Stationary Trend')&&strcmp(Model_type{3,1},'Stationary')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=avg_para(1,2:np+1)*Covariate'+avg_para(1,np+2);
                x_mu=repmat(avg_para(1,end),1,length(AMS));
            elseif strcmp(Model_type{2,1},'Non-Stationary Trend')&&strcmp(Model_type{3,1},'Non-Stationary Trend')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=avg_para(1,2:np+1)*Covariate'+avg_para(1,np+2);
                x_mu=avg_para(1,np+3:end-1)*Covariate'+avg_para(1,end);
            end
            param=[x_k;x_sigma;x_mu];
            non_xx=log(gevpdf(AMS,param(1,:),exp(param(2,:)),param(3,:)));
            L=sum(non_xx(1,:),2);
     case 'gamma'
          if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=repmat(avg_para(1,2),1,length(AMS));
          elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=avg_para(1,2:end-1)*Covariate'+avg_para(1,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
                x_k=avg_para(1,1:n)*Covariate'+avg_para(1,n+1);
                x_sigma=repmat(avg_para(1,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
               x_k=avg_para(1,1:np)*Covariate'+avg_para(1,np+1);
               x_sigma=avg_para(1,np+2:end-1)*Covariate'+avg_para(1,end);
            end
            param=[x_k;x_sigma];
            non_xx=log(gampdf(AMS,exp(param(1,:)),exp(param(2,:))));
            L=sum(non_xx(1,:),2);

     case 'gp'
            if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=repmat(avg_para(1,2),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=avg_para(1,2:end-1)*Covariate'+avg_para(1,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
                x_k=avg_para(1,1:n)*Covariate'+avg_para(1,n+1);
                x_sigma=repmat(avg_para(1,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
               x_k=avg_para(1,1:np)*Covariate'+avg_para(1,np+1);
               x_sigma=avg_para(1,np+2:end-1)*Covariate'+avg_para(1,end);
            end
            param=[x_k;x_sigma];
            non_xx=log(gppdf(AMS,(param(1,:)),exp(param(2,:))));
            L=sum(non_xx(1,:),2);

     case 'lognormal'
             if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=repmat(avg_para(1,2),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=avg_para(1,2:end-1)*Covariate'+avg_para(1,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
                x_k=avg_para(1,1:n)*Covariate'+avg_para(1,n+1);
                x_sigma=repmat(avg_para(1,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                x_k=avg_para(1,1:np)*Covariate'+avg_para(1,np+1);
                x_sigma=avg_para(1,np+2:end-1)*Covariate'+avg_para(1,end);
            end
            param=[x_k;x_sigma];
            non_xx=log(lognpdf(AMS,(param(1,:)),exp(param(2,:))));
            L=sum(non_xx(1,:),2);

        case 'normal'
             if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=repmat(avg_para(1,2),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=avg_para(1,2:end-1)*Covariate'+avg_para(1,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
                x_k=avg_para(1,1:n)*Covariate'+avg_para(1,n+1);
                x_sigma=repmat(avg_para(1,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                x_k=avg_para(1,1:np)*Covariate'+avg_para(1,np+1);
                x_sigma=avg_para(1,np+2:end-1)*Covariate'+avg_para(1,end);
            end
            param=[x_k;x_sigma];
            non_xx=log(normpdf(AMS,(param(1,:)),exp(param(2,:))));
            L=sum(non_xx(1,:),2);

     case 'ev'
             if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=repmat(avg_para(1,2),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=avg_para(1,2:end-1)*Covariate'+avg_para(1,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
                x_k=avg_para(1,1:n)*Covariate'+avg_para(1,n+1);
                x_sigma=repmat(avg_para(1,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                x_k=avg_para(1,1:np)*Covariate'+avg_para(1,np+1);
                x_sigma=avg_para(1,np+2:end-1)*Covariate'+avg_para(1,end);
            end
            param=[x_k;x_sigma];
            non_xx=log(evpdf(AMS,(param(1,:)),exp(param(2,:))));
            L=sum(non_xx(1,:),2);

     case 'weibull'
          if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=repmat(avg_para(1,2),1,length(AMS));
          elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                x_k=repmat(avg_para(1,1),1,length(AMS));
                x_sigma=avg_para(1,2:end-1)*Covariate'+avg_para(1,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
                x_k=avg_para(1,1:n)*Covariate'+avg_para(1,n+1);
                x_sigma=repmat(avg_para(1,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
               x_k=avg_para(1,1:np)*Covariate'+avg_para(1,np+1);
               x_sigma=avg_para(1,np+2:end-1)*Covariate'+avg_para(1,end);
            end
            param=[x_k;x_sigma];
            non_xx=log(wblpdf(AMS,exp(param(1,:)),exp(param(2,:))));
            L=sum(non_xx(1,:),2);

    case 'logistic'

            if strcmp(Model_type{1,1},'Stationary') && strcmp(Model_type{2,1},'Stationary')
                mu      = repmat(avg_para(1,1), 1, length(AMS));
                log_s   = repmat(avg_para(1,2), 1, length(AMS));
        
            elseif strcmp(Model_type{1,1},'Stationary') && strcmp(Model_type{2,1},'Non-Stationary Trend')
                mu      = repmat(avg_para(1,1), 1, length(AMS));
                log_s   = avg_para(1,2:end-1)*Covariate' + avg_para(1,end);
        
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend') && strcmp(Model_type{2,1},'Stationary')
                mu      = avg_para(1,1:n)*Covariate' + avg_para(1,n+1);
                log_s   = repmat(avg_para(1,end), 1, length(AMS));
        
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend') && strcmp(Model_type{2,1},'Non-Stationary Trend')
                mu      = avg_para(1,1:np)*Covariate' + avg_para(1,np+1);
                log_s   = avg_para(1,np+2:end-1)*Covariate' + avg_para(1,end);
            end
        
            % Parameter matrix (same structure)
            param = [mu; log_s];
        
            % Convert log-scale to scale
            s = exp(param(2,:));
        
            % Logistic log-likelihood
            non_xx = -log(s) ...
                     - (AMS - param(1,:))./s ...
                     - 2*log1p(exp(-(AMS - param(1,:))./s));
        
            L = sum(non_xx(1,:), 2);

 end

%% computing P

switch best_fit 
     case 'gev'
            if strcmp(Model_type{2,1},'Stationary')&& strcmp(Model_type{3,1},'Stationary')
                x_k1=repmat(para(:,1),1,length(AMS));
                x_sigma1=repmat(para(:,2),1,length(AMS));
                x_mu1=repmat(para(:,3),1,length(AMS));
            elseif strcmp(Model_type{2,1},'Stationary')&&strcmp(Model_type{3,1},'Non-Stationary Trend')
                x_k1=repmat(para(:,1),1,length(AMS));
                x_sigma1=repmat(para(:,2),1,length(AMS));
                x_mu1=para(:,3:end-1)*Covariate'+para(:,end);
            elseif strcmp(Model_type{2,1},'Non-Stationary Trend')&&strcmp(Model_type{3,1},'Stationary')
                x_k1=repmat(para(:,1),1,length(AMS));
                x_sigma1=para(:,2:np+1)*Covariate'+para(:,np+2);
                x_mu1=repmat(para(:,end),1,length(AMS));
            else
                x_k1=repmat(para(:,1),1,length(AMS));
                x_sigma1=para(:,2:np+1)*Covariate'+para(:,np+2);
                x_mu1=para(:,np+3:end-1)*Covariate'+para(:,end);
            end
            AMS1=repmat(AMS,length(para),1);
            non_xx2=log(gevpdf(AMS1,x_k1,exp(x_sigma1),x_mu1));
            sam_L=sum(non_xx2,2);
            sam_L_avg=mean(sam_L);
            P=2*(L-sam_L_avg);

     case 'gamma'
            if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
                x_k1=repmat(para(:,1),1,length(AMS));
                x_sigma1=repmat(para(:,2),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                 x_k1=repmat(para(:,1),1,length(AMS));
                 x_sigma1=para(:,2:end-1)*Covariate'+para(:,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
               x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=repmat(para(:,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
               x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=para(:,np+2:end-1)*Covariate'+para(:,end);
            end
            AMS1=repmat(AMS,length(para),1);
            non_xx2=log(gampdf(AMS1,exp(x_k1),exp(x_sigma1)));
            sam_L=sum(non_xx2,2);
            sam_L_avg=mean(sam_L);
            P=2*(L-sam_L_avg);

     case 'gp'
           if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
               x_k1=repmat(para(:,1),1,length(AMS));
               x_sigma1=repmat(para(:,2),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
               x_k1=repmat(para(:,1),1,length(AMS));
               x_sigma1=para(:,2:end-1)*Covariate'+para(:,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
               x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=repmat(para(:,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
               x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=para(:,np+2:end-1)*Covariate'+para(:,end);
            end
            AMS1=repmat(AMS,length(para),1);
            non_xx2=log(gppdf(AMS1,x_k1,exp(x_sigma1)));
            sam_L=sum(non_xx2,2);
            sam_L_avg=mean(sam_L);
            P=2*(L-sam_L_avg);

     case 'lognormal'
             if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
               x_k1=repmat(para(:,1),1,length(AMS));
               x_sigma1=repmat(para(:,2),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
               x_k1=repmat(para(:,1),1,length(AMS));
               x_sigma1=para(:,2:end-1)*Covariate'+para(:,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
               x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=repmat(para(:,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
              x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=para(:,np+2:end-1)*Covariate'+para(:,end);
            end
            AMS1=repmat(AMS,length(para),1);
            non_xx2=log(lognpdf(AMS1,x_k1,exp(x_sigma1)));
            sam_L=sum(non_xx2,2);
            sam_L_avg=mean(sam_L);
            P=2*(L-sam_L_avg);

     case 'normal'
             if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
               x_k1=repmat(para(:,1),1,length(AMS));
               x_sigma1=repmat(para(:,2),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
               x_k1=repmat(para(:,1),1,length(AMS));
               x_sigma1=para(:,2:end-1)*Covariate'+para(:,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
               x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=repmat(para(:,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
              x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=para(:,np+2:end-1)*Covariate'+para(:,end);
            end
            AMS1=repmat(AMS,length(para),1);
            non_xx2=log(normpdf(AMS1,x_k1,exp(x_sigma1)));
            sam_L=sum(non_xx2,2);
            sam_L_avg=mean(sam_L);
            P=2*(L-sam_L_avg);

        case 'ev'
             if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
               x_k1=repmat(para(:,1),1,length(AMS));
               x_sigma1=repmat(para(:,2),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
               x_k1=repmat(para(:,1),1,length(AMS));
               x_sigma1=para(:,2:end-1)*Covariate'+para(:,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
               x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=repmat(para(:,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
              x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=para(:,np+2:end-1)*Covariate'+para(:,end);
            end
            AMS1=repmat(AMS,length(para),1);
            non_xx2=log(evpdf(AMS1,x_k1,exp(x_sigma1)));
            sam_L=sum(non_xx2,2);
            sam_L_avg=mean(sam_L);
            P=2*(L-sam_L_avg);

    case 'weibull'
            if strcmp(Model_type{1,1},'Stationary')&& strcmp(Model_type{2,1},'Stationary')
                x_k1=repmat(para(:,1),1,length(AMS));
                x_sigma1=repmat(para(:,2),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Stationary')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
                 x_k1=repmat(para(:,1),1,length(AMS));
                 x_sigma1=para(:,2:end-1)*Covariate'+para(:,end);
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Stationary')
               x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=repmat(para(:,end),1,length(AMS));
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend')&&strcmp(Model_type{2,1},'Non-Stationary Trend')
               x_k1=para(:,1:np)*Covariate'+para(:,np+1);
               x_sigma1=para(:,np+2:end-1)*Covariate'+para(:,end);
            end
            AMS1=repmat(AMS,length(para),1);
            non_xx2=log(wblpdf(AMS1,exp(x_k1),exp(x_sigma1)));
            sam_L=sum(non_xx2,2);
            sam_L_avg=mean(sam_L);
            P=2*(L-sam_L_avg);

     case 'logistic'

            if strcmp(Model_type{1,1},'Stationary') && strcmp(Model_type{2,1},'Stationary')
                mu1      = repmat(para(:,1),1,length(AMS));
                log_s1   = repmat(para(:,2),1,length(AMS));
        
            elseif strcmp(Model_type{1,1},'Stationary') && strcmp(Model_type{2,1},'Non-Stationary Trend')
                mu1      = repmat(para(:,1),1,length(AMS));
                log_s1   = para(:,2:end-1)*Covariate' + para(:,end);
        
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend') && strcmp(Model_type{2,1},'Stationary')
                mu1      = para(:,1:np)*Covariate' + para(:,np+1);
                log_s1   = repmat(para(:,end),1,length(AMS));
        
            elseif strcmp(Model_type{1,1},'Non-Stationary Trend') && strcmp(Model_type{2,1},'Non-Stationary Trend')
                mu1      = para(:,1:np)*Covariate' + para(:,np+1);
                log_s1   = para(:,np+2:end-1)*Covariate' + para(:,end);
            end
        
            % Replicate data across posterior samples
            AMS1 = repmat(AMS, length(para), 1);
        
            % Convert log-scale to scale
            s1 = exp(log_s1);
        
            % Logistic log-likelihood (numerically stable)
            non_xx2 = -log(s1) ...
                      - (AMS1 - mu1)./s1 ...
                      - 2*log1p(exp(-(AMS1 - mu1)./s1));
        
            sam_L = sum(non_xx2, 2);
            sam_L_avg = mean(sam_L);
        
            P = 2*(L - sam_L_avg);


end
DIC1=-2*(L-P);
end