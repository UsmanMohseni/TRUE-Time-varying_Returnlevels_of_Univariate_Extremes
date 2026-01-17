function [var,Model_type,MW_par,Z1]=initial_var_Updated(Data,Covariate,best_fit,MW)
%%
[~,q]=size(Covariate);
switch best_fit
    case 'gev'
       var(1:q+1,1:6)=NaN;
        [MW_par]=Par_series_TSW('gev',Data,MW);
        MW_par(2,:)=log(MW_par(2,:)); % log is taken to ensure the positive values of scale parameter in GEV  
        Model_type={'-','-';'-','-';'-','-'};
        Z1(1:3,1:3)=NaN;
        for i=1:1:3
            [Z1(i,1), Z1(i,2), Z1(i,3),Model_type{i,1}]=M_K_test(MW_par(i,:));
           % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
        end
        Z1(1,end)=0;
        Model_type{1,1}='Stationary';
        %col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
       var(1,[1,3,5])=0;
        k=1;
        for i=1:1:3
           if Z1(i,3)~=0
              var(2:end,k)=0;  
           end
           k=k+2;
        end   


    case 'gp'
       var(1:q+1,1:4)=NaN;
        [MW_par]=Par_series_TSW('gp',Data,MW);
        MW_par(2,:)=log(MW_par(2,:)); % log is taken to ensure the positive values of scale parameter in GP
        Model_type={'-','-';'-','-'};
         Z1(1:2,1:3)=NaN;
        for i=1:1:2
            [Z1(i,1), Z1(i,2), Z1(i,3),Model_type{i,1}]=M_K_test(MW_par(i,:));
           % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
        end
        Z1(1,end)=0;
        Model_type{1,1}='Stationary';
        %col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
       var(1,[1,3])=0;
        k=1;
        for i=1:1:2
           if Z1(i,3)~=0
              var(2:end,k)=0;  
           end
           k=k+2;
        end        
    
    
    case 'gamma'
       var(1:q+1,1:4)=NaN;
        [MW_par]=Par_series_TSW('gam',Data,MW);
        MW_par(:,:)=log(MW_par(:,:)); % log is taken to ensure the positive values of a and b parameter in Gamma
        Model_type={'-','-';'-','-'};
        Z1(1:2,1:3)=NaN;
        for i=1:1:2
            [Z1(i,1), Z1(i,2), Z1(i,3),Model_type{i,1}]=M_K_test(MW_par(i,:));
           % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
        end
        
        %col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
       var(1,[1,3])=0;
        k=1;
        for i=1:1:2
           if Z1(i,3)~=0
              var(2:end,k)=0;  
           end
           k=k+2;
        end    


    case 'lognormal'
       var(1:q+1,1:4)=NaN;
        [MW_par]=Par_series_TSW('logn',Data,MW);
        MW_par(2,:)=log(MW_par(2,:)); % log is taken to ensure the positive values of sd parameter in LogNormal 
           Model_type={'-','-';'-','-'};
           Z1(1:2,1:3)=NaN;
        for i=1:1:2
            [Z1(i,1), Z1(i,2), Z1(i,3),Model_type{i,1}]=M_K_test(MW_par(i,:));
           % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
        end
        
        %col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
       var(1,[1,3])=0;
        k=1;
        for i=1:1:2
           if Z1(i,3)~=0
              var(2:end,k)=0;  
           end
           k=k+2;
        end  


    case 'weibull'
       var(1:q+1,1:4)=NaN;
        [MW_par]=Par_series_TSW('wbl',Data,MW);
        MW_par(:,:)=log(MW_par(:,:)); % log is taken to ensure the positive values of a and b parameter in Gamma
        Model_type={'-','-';'-','-'};
        Z1(1:2,1:3)=NaN;
        for i=1:1:2
            [Z1(i,1), Z1(i,2), Z1(i,3),Model_type{i,1}]=M_K_test(MW_par(i,:));
           % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
        end
        
        %col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
       var(1,[1,3])=0;
        k=1;
        for i=1:1:2
           if Z1(i,3)~=0
              var(2:end,k)=0;  
           end
           k=k+2;
        end 


    case 'normal'
       var(1:q+1,1:4)=NaN;
        [MW_par]=Par_series_TSW('norm',Data,MW);
        MW_par(2,:)=log(MW_par(2,:)); % log is taken to ensure the positive values of sd parameter in LogNormal 
           Model_type={'-','-';'-','-'};
           Z1(1:2,1:3)=NaN;
        for i=1:1:2
            [Z1(i,1), Z1(i,2), Z1(i,3),Model_type{i,1}]=M_K_test(MW_par(i,:));
           % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
        end
        
        %col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
       var(1,[1,3])=0;
        k=1;
        for i=1:1:2
           if Z1(i,3)~=0
              var(2:end,k)=0;  
           end
           k=k+2;
        end 


    case 'ev'
       var(1:q+1,1:4)=NaN;
        [MW_par]=Par_series_TSW('ev',Data,MW);
        MW_par(2,:)=log(MW_par(2,:)); % log is taken to ensure the positive values of sd parameter in LogNormal 
           Model_type={'-','-';'-','-'};
           Z1(1:2,1:3)=NaN;
        for i=1:1:2
            [Z1(i,1), Z1(i,2), Z1(i,3),Model_type{i,1}]=M_K_test(MW_par(i,:));
           % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
        end
        
        %col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
       var(1,[1,3])=0;
        k=1;
        for i=1:1:2
           if Z1(i,3)~=0
              var(2:end,k)=0;  
           end
           k=k+2;
        end  


        case 'logistic'
       var(1:q+1,1:4)=NaN;
        [MW_par]=Par_series_TSW('logistic',Data,MW);
        MW_par(2,:)=log(MW_par(2,:)); % log is taken to ensure the positive values of sd parameter in LogNormal 
           Model_type={'-','-';'-','-'};
           Z1(1:2,1:3)=NaN;
        for i=1:1:2
            [Z1(i,1), Z1(i,2), Z1(i,3),Model_type{i,1}]=M_K_test(MW_par(i,:));
           % col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
        end
        
        %col:1=Z value;Col:2=S value;Col:3=Trend(+1=+ve,-1=-ve,0=No trend)
       var(1,[1,3])=0;
        k=1;
        for i=1:1:2
           if Z1(i,3)~=0
              var(2:end,k)=0;  
           end
           k=k+2;
        end  

end
end
