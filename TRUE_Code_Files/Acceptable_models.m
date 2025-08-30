%% Multi-Variate Analysis CODE:2
clc
clear
close all
%% MAIN CODE
% Loading all files (xlsx format pre-cleaned data)
% Note: The file columns are as follows: Col1: Year,
% Col2=Flow,col3=Duration, col4: Volume

File_Path_AMS="E:\Mayank\2_Mahanadi\1_kantamal\1_Kantamal_Flood_Characterstics.xlsx";
File_path_Covariates="E:\Mayank\2_Mahanadi\1_kantamal\Covariate_EMD";
fileLocation="E:\Mayank\2_Mahanadi\1_kantamal\Results\Figures_EMD";
load ("E:\Mayank\2_Mahanadi\1_kantamal\Results\EMD_100_evaluations.mat");


list_Covariates=dir(fullfile(File_path_Covariates,'*.xlsx'));
[~, ~, AMS1] = xlsread(File_Path_AMS);
AMS=cell2mat(AMS1(2:end,[1,2,3,4]));


% Initialize a structure to store the data
dataStruct = struct();

% Loop through each XLSX file, load it, and store it with the file name (without extension) as the field name
for i = 1:numel(list_Covariates)
    fileName = list_Covariates(i).name;
    [~, name, ~] = fileparts(fileName);  % Remove the .xlsx extension
    filePath = fullfile(File_path_Covariates, fileName);
    dataStruct.(name) = xlsread(filePath); % Load XLSX data and store it with the modified file name as the field name
end

% Load Covariate data into separate matrices with corresponding variable names
for i = 1:numel(list_Covariates)
    fileName = list_Covariates(i).name;
    [~, name, ~] = fileparts(fileName);  % Remove the .xlsx extension
    eval([name ' = dataStruct.' name ';']);
end

% Combining all data into a master matrix
% First Column : Years
% Second Column: AMS
% rest columns : Covariates
common_years = unique(AMS(:, 1)); % Initialize with AMS years

% Initialize a cell array to store data with common years from all matrices
common_data = cell(1, numel(list_Covariates) + 3);






% Iterate through each matrix in dataStruct
for i = 1:numel(list_Covariates)
    fileName = list_Covariates(i).name;
    [~, name, ~] = fileparts(fileName);  % Remove the .xlsx extension
    current_matrix = dataStruct.(name);

    % Find common years with the current matrix
    common_years = intersect(common_years, current_matrix(:, 1));
end

% Extract and store data for the common years
common_data{1} = AMS(ismember(AMS(:, 1), common_years), 2:end); % Store AMS data

% Iterate through each matrix in dataStruct again
for i = 1:numel(list_Covariates)
    fileName = list_Covariates(i).name;
    [~, name, ~] = fileparts(fileName);  % Remove the .xlsx extension
    current_matrix = dataStruct.(name);

    % Extract data for the common years
    common_data{i + 1} = current_matrix(ismember(current_matrix(:, 1), common_years), 2);
end

common_data=cell2mat(common_data);
[m,n]=size(common_data);
Time=(1:size(common_years,1))';
DATA=nan(m,n+2);
DATA(:,1)=common_years;
DATA(:,2:end-1)=common_data;
DATA(:,end)=Time;
%% Seperating all data
Years=DATA(:,1);
PF=DATA(:,2);
PD=DATA(:,3);
PV=DATA(:,4);
Covariate=DATA(:,5:end-1);
Range=DATA(:,1);
Change_point=NaN;

%% Standardizing the DATA
[o w]=size(DATA);
DATA_Standardized(1:o,1:w)=NaN;
DATA_Standardized(:,[1,end])=DATA(:,[1,end]);

mean_PF=mean(DATA(:,2));
sd_PF=sqrt(var(DATA(:,2)));

mean_PD=mean(DATA(:,3));
sd_PD=sqrt(var(DATA(:,3)));

mean_PV=mean(DATA(:,4));
sd_PV=sqrt(var(DATA(:,4)));

DATA_Standardized(:,2)=(DATA(:,2)-mean_PF)./sd_PF;
DATA_Standardized(:,3)=(DATA(:,3)-mean_PD)./sd_PD;
DATA_Standardized(:,4)=(DATA(:,4)-mean_PV)./sd_PV;

for i=5:w-1
    avg=mean(DATA(:,i));
    SD=sqrt(var(DATA(:,i)));
     DATA_Standardized(:,i)=(DATA(:,i)-avg)./SD;%
end

%% Seperating all data
Years=DATA_Standardized(:,1);
% PF=DATA_Standardized(:,2);
% PD=DATA_Standardized(:,3);
% PV=DATA_Standardized(:,4);
Covariate=DATA_Standardized(:,5:end-1);
Range=DATA_Standardized(:,1);


[combinations]=combination_generator(size(Covariate,2));


%% best models
PF_DIC=sortrows([[1:size(combinations,1)]' combinations DIC_PF1_average],size([combinations DIC_PF1_average],2)+1);
PF_DIC(:,end+1)=PF_DIC(:,end)-PF_DIC(1,end);



PD_DIC=sortrows([[1:size(combinations,1)]' combinations DIC_PD1_average],size([combinations DIC_PD1_average],2)+1);
PD_DIC(:,end+1)=PD_DIC(:,end)-PD_DIC(1,end);

PV_DIC=sortrows([[1:size(combinations,1)]' combinations DIC_PV1_average],size([combinations DIC_PV1_average],2)+1);
PV_DIC(:,end+1)=PV_DIC(:,end)-PV_DIC(1,end);



[m,n]=size(combinations);
for i=1:1: size(combinations,1)
    
    if(PF_DIC(i,end)<=2)
       
        Best_Covariate_PF=PF_DIC(i,2:n+1);
        Best_Covariate_PF=Best_Covariate_PF(~isnan(Best_Covariate_PF));
        Best_Model_PF=PF_DIC(i,1);
        
        %% Variables parameter series
        [int1_var,Model_type1,int1_MW_par,Z_PF]=initial_var(PF,Covariate(:,Best_Covariate_PF),best_fit_PF,20);
        % [int2_var,Model_type2,int2_MW_par,Z_PD]=initial_var(PD,Covariate(:,Best_Covariate_PD),best_fit_PD,20);
        % [int3_var,Model_type3,int3_MW_par,Z_PV]=initial_var(PV,Covariate(:,Best_Covariate_PV),best_fit_PV,20);

        int1_var=Variable_parameters1(int1_var,Model_type1,Change_point,PF,int1_MW_par,Covariate(:,Best_Covariate_PF),best_fit_PF);
        % int2_var=Variable_parameters1(int2_var,Model_type2,Change_point,PD,int2_MW_par,Covariate(:,Best_Covariate_PD),best_fit_PD);
        % int3_var=Variable_parameters1(int3_var,Model_type3,Change_point,PV,int3_MW_par,Covariate(:,Best_Covariate_PV),best_fit_PV);
        % % 
        % Exponential has been taken for these coefficients
        para_PF=PARA_FINAL_PF1(:,:,Best_Model_PF);
        % para_PD=PARA_FINAL_PD1(:,:,Best_Model_PD);
        % para_PV=PARA_FINAL_PV1(:,:,Best_Model_PV);

        count_PF= size(int1_var,1)*size(int1_var,2)- sum(isnan(int1_var), 'all');
        % count_PD= size(int2_var,1)*size(int2_var,2)- sum(isnan(int2_var), 'all');
        % count_PV= size(int3_var,1)*size(int3_var,2)- sum(isnan(int3_var), 'all');

        para_PF=para_PF(:,1:count_PF);
        % para_PD=para_PD(:,1:count_PD);
        % para_PV=para_PV(:,1:count_PV);

        %% Univariate Marginals for Flood Characterstics

        [PDF_PF(:,i), CDF_PF(:,i),param_PF,params_PF,~,~]=marginal_calculator_acceptable_models(PF,para_PF,Model_type1,int1_var,Covariate(:,Best_Covariate_PF),best_fit_PF);
        % [PDF_PD, CDF_PD,para_PD,params_PD]=marginal_calculator(PD,para_PD,Model_type2,int2_var,Covariate(:,Best_Covariate_PD),best_fit_PD);
        % [PDF_PV, CDF_PV,para_PV,params_PV]=marginal_calculator(PV,para_PV,Model_type3,int3_var,Covariate(:,Best_Covariate_PV),best_fit_PV);
        %% Credible Interval
        
        if(i==1)
                    [~,~,~,~,PDF_PF_CI,CDF_PF_CI]=marginal_calculator_acceptable_models(PF,para_PF,Model_type1,int1_var,Covariate(:,Best_Covariate_PF),best_fit_PF);
        end
       
    end

end

for i=1:1: size(combinations,1)
    
    if(PD_DIC(i,end)<=2)
       
        Best_Covariate_PD=PD_DIC(i,2:n+1);
        Best_Covariate_PD=Best_Covariate_PD(~isnan(Best_Covariate_PD));
        Best_Model_PD=PD_DIC(i,1);
        
        %% Variables parameter series
        [int2_var,Model_type2,int2_MW_par,Z_PD]=initial_var(PD,Covariate(:,Best_Covariate_PD),best_fit_PD,20);

        int2_var=Variable_parameters1(int2_var,Model_type2,Change_point,PD,int2_MW_par,Covariate(:,Best_Covariate_PD),best_fit_PD);
        para_PD=PARA_FINAL_PD1(:,:,Best_Model_PD);

        
        count_PD= size(int2_var,1)*size(int2_var,2)- sum(isnan(int2_var), 'all');

        para_PD=para_PD(:,1:count_PD);

        %% Univariate Marginals for Flood Characterstics

        [PDF_PD(:,i), CDF_PD(:,i),param_PD,params_PD,~,~]=marginal_calculator_acceptable_models(PD,para_PD,Model_type2,int2_var,Covariate(:,Best_Covariate_PD),best_fit_PD);
        %% Credible Interval
        
        if(i==1)
        [~,~,~,~,PDF_PD_CI,CDF_PD_CI]=marginal_calculator_acceptable_models(PD,para_PD,Model_type2,int2_var,Covariate(:,Best_Covariate_PD),best_fit_PD);
        end
       
    end

end

for i=1:1: size(combinations,1)
    
    if(PV_DIC(i,end)<=2)
       
        Best_Covariate_PV=PV_DIC(i,2:n+1);
        Best_Covariate_PV=Best_Covariate_PV(~isnan(Best_Covariate_PV));
        Best_Model_PV=PV_DIC(i,1);
        
        %% Variables parameter series
        [int3_var,Model_type3,int3_MW_par,Z_PV]=initial_var(PV,Covariate(:,Best_Covariate_PV),best_fit_PV,20);

       int3_var=Variable_parameters1(int3_var,Model_type3,Change_point,PV,int3_MW_par,Covariate(:,Best_Covariate_PV),best_fit_PV);
        
        para_PV=PARA_FINAL_PV1(:,:,Best_Model_PV);

        
        count_PV= size(int3_var,1)*size(int3_var,2)- sum(isnan(int3_var), 'all');

        
        para_PV=para_PV(:,1:count_PV);

        %% Univariate Marginals for Flood Characterstics

     [PDF_PV(:,i), CDF_PV(:,i),param_PV,params_PV,~,~]=marginal_calculator_acceptable_models(PV,para_PV,Model_type3,int3_var,Covariate(:,Best_Covariate_PV),best_fit_PV);
        %% Credible Interval
        
        if(i==1)
                    [~,~,~,~,PDF_PV_CI,CDF_PV_CI]=marginal_calculator_acceptable_models(PV,para_PV,Model_type3,int3_var,Covariate(:,Best_Covariate_PV),best_fit_PV);
        end
       
    end

end
% 
% CDF_PF_CI=sort(CDF_PF_CI)
% CDF_PD_CI=sort(CDF_PD_CI)
% CDF_PV_CI=sort(CDF_PV_CI)
%%
figure(1)
set(gcf, 'Units', 'inches', 'Position', [1, 1, 6, 2],'color','w');

tiledlayout(1,3,"TileSpacing","loose","Padding","tight")

% Define consistent colors
mainColor = [0.1 0.4 0.8];  % Rich Blue for Best Model
confColor = [0.3 0.7 0.4];  % Pastel Green for 95% CI
modelColor = [0.5 0.5 0.5]; % Soft Gray for Acceptable Model
fillAlpha = 0.5;            % Transparency level for fills

x = linspace(1,76,76)';

%% Plot 1: CDF of PF
nexttile;
fill([x' fliplr(x')], [sort(CDF_PF_CI(:,1))' fliplr(sort(CDF_PF_CI(:,2))')], ...
     confColor, 'FaceAlpha', fillAlpha, 'EdgeColor', 'none');
hold on;
plot(sort(CDF_PF(:,1)), 'Color', mainColor, 'LineWidth', 2.5);
plot(sort(CDF_PF(:,2:end)), 'Color', modelColor, 'LineStyle', '--', 'LineWidth', 1);
title('Peak Flow', 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('Data Index', 'FontName', 'Times New Roman');
ylabel('CDF', 'FontName', 'Times New Roman');
legend('95% CI', 'Best Model', 'Acceptable Model', 'Location', 'best', 'FontName', 'Times New Roman','Fontweight','Bold');
grid on;

%% Plot 2: CDF of PD
nexttile;
fill([x' fliplr(x')], [sort(CDF_PD_CI(:,1))' fliplr(sort(CDF_PD_CI(:,2))')], ...
     confColor, 'FaceAlpha', fillAlpha, 'EdgeColor', 'none');
hold on;
plot(sort(CDF_PD(:,1)), 'Color', mainColor, 'LineWidth', 2.5);
plot(sort(CDF_PD(:,2:end)), 'Color', modelColor, 'LineStyle', '--', 'LineWidth', 1);
title('Duration', 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('Data Index', 'FontName', 'Times New Roman');
% ylabel('CDF', 'FontName', 'Times New Roman');
legend('95% CI', 'Best Model', 'Acceptable Model', 'Location', 'best', 'FontName', 'Times New Roman','Fontweight','Bold');
grid on;

%% Plot 3: CDF of PV
nexttile;
fill([x' fliplr(x')], [sort(CDF_PV_CI(:,1))' fliplr(sort(CDF_PV_CI(:,2))')], ...
     confColor, 'FaceAlpha', fillAlpha, 'EdgeColor', 'none');
hold on;
plot(sort(CDF_PV(:,1)), 'Color', mainColor, 'LineWidth', 2.5);
plot(sort(CDF_PV(:,2:end)), 'Color', modelColor, 'LineStyle', '--', 'LineWidth', 1);
title('Volume', 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('Data Index', 'FontName', 'Times New Roman');
% ylabel('CDF', 'FontName', 'Times New Roman');
legend('95% CI', 'Best Model', 'Acceptable Model', 'Location', 'best' ,'FontName', 'Times New Roman','Fontweight','Bold');
grid on;

fullPath = fullfile(fileLocation, sprintf('%s_%d', 'Acceptable Models',23)); 
    % Save in PNG format (high resolution)
    print(fullPath, '-dpng', '-r300');
    % Save in TIFF format
    print(fullPath, '-dtiff', '-r300'); 
    % Save in FIG format (MATLAB figure)
    savefig(fullPath); 
    % Save in EPS format (for vector graphics)
    hgexport(gcf, fullPath); 
    % Save in EMF format (for compatibility with MS Office)
    print(fullPath, '-dmeta');  

