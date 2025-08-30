%% Function Variable Parameters
function [int_var]= Variable_parameters(int_var,Model_type,change_point,AMS,par,Covariate)
Years=size(AMS,1);
Covariate1=Covariate(11:Years-9,:);

[p,~]=size(Model_type);
for i=1:1:p
if(strcmp(Model_type{i,1},'Stationary'))
    p=gevfit(AMS);
    int_var(1,2*i-1)=p(:,i);
elseif(strcmp(Model_type{i,1},'Non-Stationary Trend'))
   Model=fitlm(Covariate1,par(i,:)');
   coefficients = Model.Coefficients.Estimate;
    B=[coefficients(2:end,1);coefficients(1,1)];
   coefficients=B;
   % int_var(:,2*i-1)=flipud(coefficients);
 int_var(:,2*i-1)=coefficients;
elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step Trend'))
    p1=polyfit([11:10+change_point],par(i,1:change_point),1);
    int_var(:,2*i-1)=p1;
    p2=polyfit([11:10+change_point],par(i,change_point+1:end),1);
    int_var(:,2*i)=p2;
elseif(strcmp(Model_type{i,1},'Step Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
    p1=polyfit([11:10+change_point],par(i,1:change_point),1);
    int_var(:,2*i-1)=p1;
    p2=gevfit(AMS(change_point+1:end,1));
    int_var(1,2*i)=p2(:,i);
elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step Trend'))
     
    p1=gevfit(AMS(1:change_point,1));
    int_var(1,2*i-1)=p1(:,i);
    p2=polyfit([11:10+change_point],par(i,change_point+1:end,1),1);
    int_var(:,2*i)=p2;
elseif(strcmp(Model_type{i,1},'Step No-Trend')&&strcmp(Model_type{i,2},'Step No-Trend'))
    p1=gevfit(AMS(1:change_point,1));
    int_var(1,2*i-1)=p1(:,i);
    p2=gevfit(AMS(change_point+1:end,1));
    int_var(1,2*i)=p2(:,i);
end
end
end
