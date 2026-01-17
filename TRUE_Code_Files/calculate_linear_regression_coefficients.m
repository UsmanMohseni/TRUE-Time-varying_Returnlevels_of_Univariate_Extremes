%% Function for calculating coefficients
function coefficients = calculate_linear_regression_coefficients(covariates, target)
    % Perform linear regression
    data = [covariates, target];
    
    % Extract variable names from the input covariates matrix
    num_covariates = size(covariates, 2);
    variable_names = cell(1, num_covariates);
    for i = 1:num_covariates
        variable_names{i} = ['Covariate', num2str(i)];
    end
    variable_names{end + 1} = 'Z';
    
    tbl = array2table(data, 'VariableNames', variable_names);
    formula = 'Z ~ ';
    
    for i = 1:num_covariates
        formula = [formula, variable_names{i}, ' + '];
    end
    formula = formula(1:end-2); % Remove the trailing '+'
    
    lm = fitlm(tbl, formula);

    % Extract and return the coefficients
    coefficients = lm.Coefficients.Estimate;
end
