function [param]=validate_params(para,Model_type,best_fit,Covariate)
   [m,n]=size(Covariate);
    switch best_fit
        case 'gev'
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
        case 'gp'
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
        case 'gamma'
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
        case 'lognormal'
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
        case 'frank'
             if(strcmp(Model_type,'Stationary'))
                teta=repmat(para(:,1),1,m);
             elseif(strcmp(Model_type,'Non-Stationary Trend'))
                teta=para(:,1:n)*Covariate' + para(:,n+1);
             end
             param=median(teta);
    end
end