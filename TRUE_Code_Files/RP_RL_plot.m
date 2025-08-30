function RP_RL_plot(scale_t,shape,mu_t,Time,sd1,mu1,no,Label)
    X=ceil(prctile(Time,95));
    scale_t=scale_t(:,X);
    shape=shape(:,X);
    mu_t=mu_t(:,X);
    para=[shape scale_t mu_t];
    z=norminv(1-0.05/2);
    CI(1:3,1:2)=NaN;
    
    for i=1:1:3
        mu(i,1)=mean(para(:,i));
        
        sd=std(para(:,i));
        CI(i,1)=mu(i,1)-z*sd;
        CI(i,2)=mu(i,1)+z*sd;
    end
z2=norminv(1-0.01/2);
    CI2(1:3,1:2)=NaN;
    
    for i=1:1:3
        mu(i,1)=mean(para(:,i));
        
        sd=std(para(:,i));
        CI2(i,1)=mu(i,1)-z2*sd;
        CI2(i,2)=mu(i,1)+z2*sd;
    end
    
    
   figure(no)

    for i=2:1:100
        RL1(i-1,:)=(gevinv(1-1/i,mu(1,1),mu(2,1),mu(3,1)))*sd1+mu1;
        RL2(i-1,:)=(gevinv(1-1/i,CI(1,1),CI(2,1),CI(3,1)))*sd1+mu1;
        RL3(i-1,:)=(gevinv(1-1/i,CI(1,2),CI(2,2),CI(3,2)))*sd1+mu1;

    end
    for i=2:1:100
        
        RL4(i-1,:)=(gevinv(1-1/i,CI2(1,1),CI2(2,1),CI2(3,1)))*sd1+mu1;
        RL5(i-1,:)=(gevinv(1-1/i,CI2(1,2),CI2(2,2),CI2(3,2)))*sd1+mu1;

    end
    
    X=[2:100]';
    
%     [k,l]=size(RL);
%     for i=1:1:k
%         RL1(i,:)=prctile(RL(i,:),50);
%         RL2(i,:)=prctile(RL(i,:),95);
%         RL3(i,:)=prctile(RL(i,:),5);
%         RL4(i,:)=prctile(RL(i,:),0);
%         RL5(i,:)=prctile(RL(i,:),100);
%     end
     fill([2:100,fliplr(2:100)],[RL4',fliplr(RL5')],'y','edgecolor','none')
     hold on
   
    plot(X,RL1,'color',[0.5 0 0],'LineWidth',4)
    hold on
    
    plot(X,RL2,'-.','color',[0 0 0],'LineWidth',2)
    hold on
   
         plot(X,RL3,'-.','color',[0 0 0],'Linewidth',2)
    hold on
   

    title('Return Period Vs Return Level')
    xlabel('Return Period (in years)')
    ylabel(Label)
     legend('Ensembles','Median','95 %','5 %','location','southoutside','Orientation','horizontal')
   
end