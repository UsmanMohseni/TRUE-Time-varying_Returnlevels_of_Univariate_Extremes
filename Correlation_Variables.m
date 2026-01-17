
clc; clear all;

% Read the Excel file
data = readtable('E:\SERB_Project\2_baramghat\Streamflow\2_cleaned_Barmanghat_cleaned\Mumbai.xlsx');

% Extract individual columns
P = data{:,2};     % Precipitation
Q = data{:,3};     % Streamflow (Discharge)
TWL = data{:,4};   % Total Water Level

% Pearson Correlation
pearson_P_Q = corr(P, Q, 'Type', 'Pearson');
pearson_Q_TWL = corr(Q, TWL, 'Type', 'Pearson');
pearson_P_TWL = corr(P, TWL, 'Type', 'Pearson');

% Spearman Correlation
spearman_P_Q = corr(P, Q, 'Type', 'Spearman');
spearman_Q_TWL = corr(Q, TWL, 'Type', 'Spearman');
spearman_P_TWL = corr(P, TWL, 'Type', 'Spearman');

% Kendall Correlation
kendall_P_Q = corr(P, Q, 'Type', 'Kendall');
kendall_Q_TWL = corr(Q, TWL, 'Type', 'Kendall');
kendall_P_TWL = corr(P, TWL, 'Type', 'Kendall');

% Display results
fprintf('\n--- Pearson Correlation ---\n');
fprintf('P & Q   : %.4f\n', pearson_P_Q);
fprintf('Q & TWL : %.4f\n', pearson_Q_TWL);
fprintf('P & TWL : %.4f\n', pearson_P_TWL);

fprintf('\n--- Spearman Correlation ---\n');
fprintf('P & Q   : %.4f\n', spearman_P_Q);
fprintf('Q & TWL : %.4f\n', spearman_Q_TWL);
fprintf('P & TWL : %.4f\n', spearman_P_TWL);

fprintf('\n--- Kendall Correlation ---\n');
fprintf('P & Q   : %.4f\n', kendall_P_Q);
fprintf('Q & TWL : %.4f\n', kendall_Q_TWL);
fprintf('P & TWL : %.4f\n', kendall_P_TWL);
