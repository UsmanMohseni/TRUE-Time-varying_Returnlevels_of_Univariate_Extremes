%% Function for MK test

function [Z, S ,C,Trend] = M_K_test (X)
    S=0;
    X=X';
    n=size(X,1);
    for i=1:1:n
        for j=i+1:1:n
            S=S+sign(X(j,1)-X(i,1));
        end
    end
    varS = (n)*(n-1)*(2*n+5)/18;

    if S>0
        Z=(S-1)/sqrt(varS);
    elseif S == 0
        Z = 0;
    elseif S<0
        Z=(S+1)/sqrt(varS);
    end

     if (abs(Z)>=1.65 && S<0)
             C=-1;
             Trend='Non-Stationary Trend';
         elseif(abs(Z)>=1.65 && S>0)
             C=1;
             Trend='Non-Stationary Trend';
     elseif(abs(Z)<1.65)
             C=0;
             Trend='Stationary';
    end
         
end