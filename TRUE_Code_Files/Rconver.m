%% Function Rconver
function Rhat= Rconver(ysi,evl,bur,cha,mix,sts,sur,dim)
% check convergence using R hat by Gelman
for Rsigm= ysi
    if strcmp(sur,'sta');
        for col= 1:dim
            pa(:,col,:)= mix(:,col,bur:evl);
            Rv_pa(col,:)= reshape(pa(:,col,:),1,cha*(evl-bur+1),1);
            var_pa(:,col)= var(reshape(pa(:,col,:),cha,(evl-bur+1)),0,2);
            i=var()
            ave_pa(:,col)= mean(var_pa(:,col),1);
            var_mix_pa(:,col)= var(Rv_pa(col,sts:end),0,2);
            Rhat(:,col)= sqrt(var_mix_pa(:,col)./ave_pa(:,col));
            if Rhat(:,col)> 1.1
                S=sprintf('Station %i, %i parameter out of Rhat bound', Rsigm, col);
                disp(S)
                disp('Maybe increase the evaluation')
            end
        end
        Rhat= Rhat;
    else
        for coll= 1:dim
            pa(:,coll,:)= mix(:,coll,bur:evl);
            Rv_pa(coll,:)= reshape(pa(:,coll,:),1,cha*(evl-bur+1),1);
            var_pa(:,coll)= var(reshape(pa(:,coll,:),cha,(evl-bur+1)),0,2);
            ave_pa(:,coll)= mean(var_pa(:,coll),1);
            var_mix_pa(:,coll)= var(Rv_pa(coll,sts:end),0,2);
            Rhat(:,coll)= sqrt(var_mix_pa(:,coll)./ave_pa(:,coll));
            if Rhat(:,coll)> 1.1
                S=sprintf('Station %i, %i parameter out of Rhat bound', Rsigm, coll);
                disp(S)
                disp('Maybe increase the evaluation')
            end
        end
        Rhat= Rhat;
    end
end
end