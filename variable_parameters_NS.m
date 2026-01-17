function [shape, scale_t, mu_t ]=variable_parameters_NS(Model_type,para,Covariate)
%This code is without change point
[m,n]=size(Covariate);
if(strcmp(Model_type(2,1),'Stationary') && strcmp(Model_type(3,1),'Stationary'))
    shape=repmat(para(:,1),1,m);
    scale_t=repmat(para(:,2),1,m);
    mu_t=repmat(para(:,3),1,m);
elseif(strcmp(Model_type(2,1),'Stationary') && strcmp(Model_type(3,1),'Non-Stationary Trend'))
    shape=repmat(para(:,1),1,m);
    scale_t = repmat(para(:,2),1,m);
    mu_t = para(:,3:end-1)*Covariate' + para(:,end);
    % mu_t = para(:,4:end)*Covariate' + para(:,3);
elseif(strcmp(Model_type(2,1),'Non-Stationary Trend') && strcmp(Model_type(3,1),'Stationary'))
    shape=repmat(para(:,1),1,m);
    scale_t= para(:,2:n+1)*Covariate' + para(:,n+2);
    mu_t = repmat(para(:,end),1,m);
elseif(strcmp(Model_type(2,1),'Non-Stationary Trend') && strcmp(Model_type(3,1),'Non-Stationary Trend'))
    shape=repmat(para(:,1),1,m);
    scale_t= para(:,2:n+1)*Covariate' + para(:,n+2);
    mu_t = para(:,n+3:end-1)*Covariate'+ para(:,end);
end
scale_t=exp(scale_t);
% mu11=median(mu_t);
% shape11=median(shape);
% scale11=median(scale_t);
% param=[shape11; scale11; mu11];
end