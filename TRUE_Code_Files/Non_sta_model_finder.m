clc
clear
close all

load ("E:\Mayank\Files\DATA\WRIS_Dataset_processed\Dataset\BaseflowSeperatedDataWithCharacterstics.mat")




%% Finding Best Distribution fit
% fit_PF=fitmethis2(PF);
% best_fit_PF=fit_PF(1).name;
% fit_PD=fitmethis2(PD);
% best_fit_PD=fit_PD(1).name;
% fit_PV=fitmethis2(PV);
% best_fit_PV=fit_PV(1).name;
% [int1_var,Model_type1,int1_MW_par]=initial_var(PF,Covariate(:,Best_Covariate_PF),best_fit_PF,20);
% [int2_var,Model_type2,int2_MW_par]=initial_var(PD,Covariate(:,Best_Covariate_PD),best_fit_PD,20);
% [int3_var,Model_type3,int3_MW_par]=initial_var(PV,Covariate(:,Best_Covariate_PV),best_fit_PV,20);

Data=BaseflowSeperatedData.FloodCharacterstics;
Model_type1 = cell(268, 3); % Initialize a 268x3 cell array
Model_type1(1:268, 1:3) = {'NA'}; % Assign 'NA' to all elements
Model_type2 = cell(268, 3); % Initialize a 268x3 cell array
Model_type2(1:268, 1:3) = {'NA'}; % Assign 'NA' to all elements
Model_type3 = cell(268, 3); % Initialize a 268x3 cell array
Model_type3(1:268, 1:3) = {'NA'}; % Assign 'NA' to all elements

warning("off")
for i=1:1:length(Data)

    Cts=Data{i};
    years=Cts(:,1);
    sizey(i,1)=length(years);
    PF=Cts(:,2);
    PD=Cts(:,3);
    PV=Cts(:,4);
     fit_PF=fitmethis2(PF);
    best_fit_PF{i,1}=fit_PF(1).name;
    fit_PD=fitmethis2(PD);
    best_fit_PD{i,1}=fit_PD(1).name;
    fit_PV=fitmethis2(PV);
    best_fit_PV{i,1}=fit_PV(1).name;
    if(length(years)>=35)
    %% Finding Best Distribution fit
   
    try

    [~,Model_type11,~]=initial_var(PF,years,best_fit_PF{i,1},20);
    if(size(Model_type11,1)==2)
        Model_type1(i,1:2)=Model_type11(:,1)';
    else
        Model_type1(i,:)=Model_type11(:,1)';
    end
    [~,Model_type22,~]=initial_var(PD,years,best_fit_PD{i,1},20);
     if(size(Model_type22,1)==2)
        Model_type2(i,1:2)=Model_type22(:,1)';
    else
        Model_type2(i,:)=Model_type22(:,1)';
    end
    [~,Model_type33,~]=initial_var(PV,years,best_fit_PV{i,1},20);
     if(size(Model_type33,1)==2)
        Model_type3(i,1:2)=Model_type33(:,1)';
    else
        Model_type3(i,:)=Model_type33(:,1)';
    end
    catch msg
        fprintf(msg.message)
        % best_fit_PF{i,1}='NA';
        %  best_fit_PD{i,1}='NA';
        %  best_fit_PV{i,1}='NA';
        % Model_type1{i,1}='NA';
        % Model_type2{i,1}='NA';
        % Model_type3{i,1}='NA';
    end
    else
         % best_fit_PF{i,1}='NA';
         % best_fit_PD{i,1}='NA';
         % best_fit_PV{i,1}='NA';
        % Model_type1{i,1}='NA';
        % Model_type2{i,1}='NA';
        % Model_type3{i,1}='NA';
        
    end
disp(i)
end