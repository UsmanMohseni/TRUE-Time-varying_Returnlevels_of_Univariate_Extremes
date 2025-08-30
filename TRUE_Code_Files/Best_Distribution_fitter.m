clc
clear
close all
% Load your data

File_Path_AMS="D:\MTech\Narmada_Basin\Flood_characterstics\Handia.xlsx";
[~, ~, AMS1] = xlsread(File_Path_AMS);
AMS=cell2mat(AMS1(2:end,[1,4]));
data=AMS(:,2);
% Define candidate distributions
distributions = {'gamma', 'exponential','gev','GeneralizedPareto'};

% Initialize arrays to store AIC and BIC values
AIC = zeros(size(distributions));
BIC = zeros(size(distributions));

% Fit each distribution to the data and compute AIC and BIC
for i = 1:length(distributions)
    dist = fitdist(data, distributions{i});
    params = dist.ParameterValues;
    n = numel(data);
    log_likelihood = sum(log(pdf(dist, data)));
    k = numel(params);
    AIC(i) = -2 * log_likelihood + 2 * k;
    BIC(i) = -2 * log_likelihood + k * log(n);
end

% Choose the distribution with the minimum AIC or BIC
[minAIC, minAICIdx] = min(AIC);
[minBIC, minBICIdx] = min(BIC);

bestFitAIC = distributions{minAICIdx};
bestFitBIC = distributions{minBICIdx};

disp(['Best fit distribution using AIC: ' bestFitAIC]);
disp(['Best fit distribution using BIC: ' bestFitBIC]);
AIC=AIC';
BIC=BIC';
