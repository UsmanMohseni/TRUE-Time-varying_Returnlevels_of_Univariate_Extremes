%-------------------------------------------------------------------------
%                           RUN TRUE_NSMCFFA
%-------------------------------------------------------------------------

%Reference Publication: Mohseni U, Gailakoti M, Vinnarasi R, 2025
% Non-Stationary Multi-Covariate Flood Frequency Analysis with TRUE: 
% A Tool for Time-varying Return levels of Univariate Extremes
% 

%TRUE: Time-varying Return level for Univariate Extremes

% TRUE MATLAB code is developed by Usman Mohseni, Research Scholar, 
% Department of Civil Engineering,
% Indian Institute of Technology Roorkee, Uttarakhand, India.

% For questions and permissions, please contact 
% Dr. Vinnarasi Rajendran: vinnarasi@ce.iitr.ac.in or
% Usman Mohseni: mohseni_ua@ce.iitr.ac.in 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disclaimer: The TRUE (Time-varying Return level for Univariate Extremes) 
% toolbox is made available freely for research and educational purposes. 
% It is distributed without any warranties or guarantees, either expressed 
% or implied, regarding the accuracy, completeness, or suitability of the 
% code, results, or visualizations. Users are fully responsible for how 
% they apply TRUE and for any consequences that may arise from its use. 
% The developers and affiliated institutions shall not be held liable for 
% any loss, damage, or costs (whether direct, indirect, incidental, or 
% consequential) connected to the use of this software. The algorithms and 
% features may be updated or modified without prior notice.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all


%% Enter the location of data files

tic
File_Path_AMS= "E:\TRUE\Data\Streamflow\4_cleaned_Handia_cleaned_Analysis.xlsx";
File_path_Covariates="E:\TRUE\Data\Covariates";
fileLocation = "E:\TRUE"; % location to store the figures
finalsave = "E:\TRUE\All_Results.mat"; % workspace save location

Total_simulations = 20; % Total number of simulations (modifiable by the user as required)


%% Load AMS data and covariates

list_Covariates = dir(fullfile(File_path_Covariates, '*.xlsx'));
[~, ~, AMS1] = xlsread(File_Path_AMS);
PF = cell2mat(AMS1(2:end, [1,2])); % First col = years, second col = PF values

% Store covariates in a structure
dataStruct = struct();
for i = 1:numel(list_Covariates)
    fileName = list_Covariates(i).name;
    [~, name, ~] = fileparts(fileName);
    filePath = fullfile(File_path_Covariates, fileName);
    dataStruct.(name) = xlsread(filePath);
end

% Merge PF and covariates for common years
common_years = PF(:,1);
common_data = cell(1, numel(list_Covariates) + 1);

common_data{1} = PF(:,2); % PF values
for i = 1:numel(list_Covariates)
    fileName = list_Covariates(i).name;
    [~, name, ~] = fileparts(fileName);
    current_matrix = dataStruct.(name);
    % Intersect years
    common_years = intersect(common_years, current_matrix(:,1));
end

% Extract for common years
common_data{1} = PF(ismember(PF(:,1), common_years), 2);
for i = 1:numel(list_Covariates)
    fileName = list_Covariates(i).name;
    [~, name, ~] = fileparts(fileName);
    current_matrix = dataStruct.(name);
    common_data{i+1} = current_matrix(ismember(current_matrix(:,1), common_years), 2);
end

% Final DATA matrix
common_data = cell2mat(common_data);
[m, n] = size(common_data);
Time = (1:size(common_years,1))';
DATA = nan(m, n+2);
DATA(:,1) = common_years;
DATA(:,2:end-1) = common_data;
DATA(:,end) = Time;

%% Separate PF and covariates

Years = DATA(:,1);
PF = DATA(:,2);
Covariate = zscore(DATA(:,3:end-1)); % Standardized covariates
Covariate1 = Covariate;

%% Model setup

numVariables = size(Covariate, 2);
allCombinations = cell(1, numVariables);
for variableIndex = 1:numVariables
    allCombinations{variableIndex} = nchoosek(1:numVariables, variableIndex);
end

PARA_FINAL_PF1(1:25005, 1:15, 1:63) = 0;

%% Find best distribution fit (exclude Gamma if needed)
fit_PF=fitmethis2(PF);
fit_PF = fit_PF(~strcmpi({fit_PF.name}, 'gamma')); % remove gamma
best_fit_PF=fit_PF(1).name;


% fit_PF = fitmethis2(PF);
% fit_PF = fit_PF(~strcmpi({fit_PF.name}, 'gamma')); % remove gamma
% % fit_PF = fit_PF(~strcmpi({fit_PF.name}, 'gp')); % remove gamma
% [~, idx] = min([fit_PF.aic]);
% best_fit_PF = fit_PF(idx).name;

%% Simulations for DIC evaluation

c = parcluster('local');
parpool(c, c.NumWorkers);
disp('Simulations are running....')

parfor simulation = 1:Total_simulations
    [DIC_PF1(:,simulation),~,PARA_FINAL_PF] = model_evaluator2(numVariables, allCombinations, Covariate1, PF, best_fit_PF);
    PARA_FINAL_PF1 = PARA_FINAL_PF1 + PARA_FINAL_PF;

    disp(['Simulation ' num2str(simulation) ' done']);
end

PARA_FINAL_PF1 = PARA_FINAL_PF1 ./ Total_simulations;
DIC_PF1_average = mean(DIC_PF1, 2);

delete(gcp('nocreate'));

%% Identify best covariates

n = size(Covariate1,2);
combinations = {};
index = 1;
for k = 1:n
    combs = nchoosek(1:n, k);
    for i = 1:size(combs,1)
        combinations{index} = combs(i,:);
        index = index + 1;
    end
end

[~, BM_PF_index] = min(DIC_PF1_average);
Best_Covariate_PF = cell2mat(combinations(BM_PF_index));
Best_Model_PF = BM_PF_index;

%% Standardize PF and covariates

[o, w] = size(DATA);
DATA_Standardized = NaN(o, w);
DATA_Standardized(:,[1,end]) = DATA(:,[1,end]);

mean_PF = mean(DATA(:,2));
sd_PF = std(DATA(:,2));
DATA_Standardized(:,2) = (DATA(:,2)-mean_PF)./sd_PF;

for i = 3:w-1
    avg = mean(DATA(:,i));
    SD = std(DATA(:,i));
    DATA_Standardized(:,i) = (DATA(:,i)-avg)./SD;
end

%% Variables parameter series

Years = DATA_Standardized(:,1);
Covariate = DATA_Standardized(:,3:end-1);
Change_point = NaN;

[int1_var,Model_type1,int1_MW_par,Z_PF] = initial_var(PF, Covariate(:,Best_Covariate_PF), best_fit_PF, 20);
int1_var = Variable_parameters1(int1_var, Model_type1, Change_point, PF, int1_MW_par, Covariate(:,Best_Covariate_PF), best_fit_PF);

para_PF = PARA_FINAL_PF1(:,:,Best_Model_PF);
count_PF = numel(int1_var) - sum(isnan(int1_var),'all');
para_PF = para_PF(:,1:count_PF);

[PDF_PF, CDF_PF, para_PF, params_PF] = marginal_calculator(PF, para_PF, Model_type1, int1_var, Covariate(:,Best_Covariate_PF), best_fit_PF);

%% Stationary modelling

cha = 5; evl = 15000; bur = 10000; sts = 1;
[sta_int1_var,sta_Model_type1,sta_int1_MW_par,~] = initial_var(PF, Covariate(:,1), best_fit_PF, 20);
sta_int1_var(2,:) = NaN;
sta_Model_type1(:,1) = {'Stationary'};
sta_int1_var = Variable_parameters1(sta_int1_var, sta_Model_type1, Change_point, PF, sta_int1_MW_par, Covariate(:,Best_Covariate_PF), best_fit_PF);

[sta_para_PF,sta_Rhat_PF,~,~] = demc_all(sta_Model_type1, sta_int1_var, PF, evl, bur, sts, cha, Covariate(:,1), NaN, best_fit_PF);
while max(sta_Rhat_PF) > 1.1 || any(isnan(sta_Rhat_PF))
    [sta_para_PF,sta_Rhat_PF,~,~] = demc_all(sta_Model_type1, sta_int1_var, PF, evl, bur, sts, cha, Covariate(:,1), NaN, best_fit_PF);
end

[sta_PDF_PF, sta_CDF_PF,~,sta_params_PF] = marginal_calculator_without_log(PF, sta_para_PF, sta_Model_type1, sta_int1_var, Covariate(:,1), best_fit_PF);

%% Save results

save(finalsave,'-v7.3');
toc


%% Plots 

figure(1);
set(gcf, 'Color', 'w', 'Units', 'centimeters', 'Position', [10, 4.5, 18.75, 9.15]); % Set size to 8x10 inches

DATA = sortrows(DATA,1);

% Use only available years (no missing year handling)
DataRange = DATA(:,1);

% Define colors for each plot
colors = [
    0.1, 0.4, 0.8; % Blue for Peak Duration
    0.3, 0.6, 0.1;  % Green for Peak Discharge 
    0.9, 0.3, 0.3   % Red for Peak Volume
];

lineWidth = 2; % Thicker line width for better visibility
gridColor = [0.9, 0.9, 0.9]; % Very Light Gray for grid
markerStyle = 'o'; % Circle markers
markerSize = 5; % Size of markers

% Set font properties for the entire figure
set(groot, 'DefaultAxesFontName', 'Times New Roman'); % Change font for axes
set(groot, 'DefaultTextFontName', 'Times New Roman'); % Change font for text

% Plot Discharge Time Series

subplot(1, 1, 1);
plot(DataRange, DATA(:, 2), 'Color', colors(1, :), 'LineWidth', lineWidth, ...
    'Marker', markerStyle, 'MarkerFaceColor', colors(1, :), 'MarkerSize', markerSize);
xlim([min(DataRange) max(DataRange)]);
xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Discharge (m^3/sec)', 'FontWeight', 'bold', 'FontSize', 12);
grid on; % Add grid
set(gca, 'GridColor', gridColor); % Set grid color
set(gca, 'FontSize', 10); % Font size for axis labels
set(gca, 'LooseInset', get(gca, 'TightInset') * 1.1); % Adjust margins

fullPath = fullfile(fileLocation, sprintf('%s_%d', 'Discharge_Time_Series', 1)); 
    % Save in PNG format (high resolution)
    print(fullPath, '-dpng', '-r400');
    % Save in TIFF format
    print(fullPath, '-dtiff', '-r400'); 
    % Save in FIG format (MATLAB figure)
    savefig(fullPath); 
    % Save in EPS format (for vector graphics)
    hgexport(gcf, fullPath); 
    % Save in EMF format (for compatibility with MS Office)
    print(fullPath, '-dmeta');  


    %% Trace Plots

Trace_Plot_Univariate(para_PF, Model_type1, int1_var, best_fit_PF,2)

fullPath = fullfile(fileLocation, sprintf('%s_%d', 'Trace_plot', 2)); 
    % Save in PNG format (high resolution)
    print(fullPath, '-dpng', '-r400');
    % Save in TIFF format
    print(fullPath, '-dtiff', '-r400'); 
    % Save in FIG format (MATLAB figure)
    savefig(fullPath); 
    % Save in EPS format (for vector graphics)
    hgexport(gcf, fullPath); 
    % Save in EMF format (for compatibility with MS Office)
    print(fullPath, '-dmeta');  


    %% Trend Plot 

int_MW_par = NaN(size(int1_MW_par,1), size(int1_MW_par,2), 1);
int_MW_par(:,:,1) = int1_MW_par;

% Call plotting function only for PF
Univariate_Parameter_Trend_Plot(int_MW_par, {best_fit_PF}, Years, 20, 3)

% Define save path
fullPath = fullfile(fileLocation, sprintf('%s_%d', 'Parameter_Trend_plot', 3)); 

% Save in multiple formats
print(fullPath, '-dpng', '-r400');    % PNG high resolution
print(fullPath, '-dtiff', '-r400');   % TIFF
savefig(fullPath);                    % FIG (MATLAB figure)
hgexport(gcf, fullPath);              % EPS (vector graphics)
print(fullPath, '-dmeta');            % EMF (MS Office compatibility)



%% Return Period Vs Return Level Plots

Stationary_Univariate_RP_RL_plot(sta_params_PF,best_fit_PF,Time,'Return Level (Discharge m^3/s)',4)

fullPath = fullfile(fileLocation, sprintf('%s_%d', 'S_RP_RL_plot_univariate',4)); 
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



%  Non-Stationary Univariate

Non_Stationary_Univariate_RP_RL_plot(params_PF,best_fit_PF,Time,'Return Level (Discharge m^3/s)',5)

fullPath = fullfile(fileLocation, sprintf('%s_%d', 'NS_RP_RL_plot_univariate',5)); 
    % Save in PNG format (high resolution)
    print(fullPath, '-dpng', '-r400');
    % Save in TIFF format
    print(fullPath, '-dtiff', '-r400'); 
    % Save in FIG format (MATLAB figure)
    savefig(fullPath); 
    % Save in EPS format (for vector graphics)
    hgexport(gcf, fullPath); 
    % Save in EMF format (for compatibility with MS Office)
    print(fullPath, '-dmeta');  
