function [PDF,CDF,para,params]=marginal_calculator(DATA,para,Model_type,int_var,Covariate,best_fit)
[m,n]=size(Covariate);
    switch best_fit
        case 'gev'
            if(strcmp(Model_type(2,1),'Stationary'))
                column=2;
            elseif(strcmp(Model_type(2,1),'Non-Stationary Trend'))
                column=2:size(int_var,1);   
            end
            % para(:,column)=log(para(:,column));

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

            params{1}=shape;
            params{2}=scale_t;
            params{3}=mu_t;
            
            PDF=gevpdf(DATA,param(1,:)',exp(param(2,:))',param(3,:)');
            CDF=gevcdf(DATA,param(1,:)',exp(param(2,:))',param(3,:)');

        case 'gp'
            if(strcmp(Model_type(2,1),'Stationary'))
                column=2;
            elseif(strcmp(Model_type(2,1),'Non-Stationary Trend'))
                column=2:size(int_var,1);   
            end
            % para(:,column)=log(para(:,column));

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
            
            params{1}=shape;
            params{2}=scale_t;
         
            PDF=gppdf(DATA,(param(1,:))',exp(param(2,:))');
            CDF=gpcdf(DATA,exp(param(1,:))',exp(param(2,:))');
          
        case 'gamma'
            % para(:,:)=log(para(:,:));
           
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
              params{1}=shape;
            params{2}=scale_t;
              PDF=gampdf(DATA,exp(param(1,:))',exp(param(2,:))');
              CDF=gamcdf(DATA,exp(param(1,:))',exp(param(2,:))');

        case 'lognormal'
            if(strcmp(Model_type(2,1),'Stationary'))
                column=2;
            elseif(strcmp(Model_type(2,1),'Non-Stationary Trend'))
                column=2:size(int_var,1);   
            end
            % para(:,column)=log(para(:,column));
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
               params{1}=shape;
            params{2}=scale_t;

              PDF=lognpdf(DATA,(param(1,:))',exp(param(2,:))');
              CDF=logncdf(DATA,(param(1,:))',exp(param(2,:))');

    end
    
end