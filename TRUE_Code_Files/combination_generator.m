function [matrix]=combination_generator(n)
% n=No. of Covariates
        % Define the set of variables
        variables =1:n;
        
        % Initialize an empty cell to store combinations
        all_combinations = {};
        
        % Loop through all possible combination lengths (from 1 to 5)
        for k = 1:length(variables)
            comb = nchoosek(variables, k);  % Generate k-combinations
            all_combinations = [all_combinations; num2cell(comb, 2)];  % Store in cell array
        end
        
        % Convert the cell array to a matrix (optional: pad with NaNs to keep dimensions equal)
        max_len = max(cellfun(@length, all_combinations));
        matrix = cellfun(@(x) [x, NaN(1, max_len - length(x))], all_combinations, 'UniformOutput', false);
        matrix = cell2mat(matrix);
end