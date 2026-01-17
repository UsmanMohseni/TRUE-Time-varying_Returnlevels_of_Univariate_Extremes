function trace_plot(para,Model_type,int_var,No)
[~,n]=size(para);
z=1;
y=1;

figure(No)
switch best_dist 
    case 'gev'
if(strcmp(Model_type(1,1),'Stationary'))
    subplot(n,2,z)
    plot(para(:,y));
    axis tight;
    xlabel('Iteration')
    ylabel('\gamma')
    hold on
    z=z+1;
    subplot(n,2,z)
    hold on
    pd = fitdist(para(:,y), 'Normal'); 
    x = linspace(min(para(:,y)), max(para(:,y)), 100);
    y_fit = pdf(pd,x);
    histogram(para(:,y),100,'Normalization','pdf','FaceColor','b','EdgeColor','none');
    plot(x, y_fit, 'r', 'LineWidth', 2); % Plot the curve in red
    hold on
    % xlabel('\gamma')
    ylabel('Frequency')
    z=z+1;
    y=y+1;
end

if(strcmp(Model_type(1,1),'Non-Stationary Trend'))
    subplot(n,2,z)
    plot(para(:,y));
    axis tight;
    xlabel('Iteration')
    ylabel('\gamma_t_1')
    hold on
    z=z+1;
    subplot(n,2,z)
    hold on
    pd = fitdist(para(:,y), 'Normal'); 
    x = linspace(min(para(:,y)), max(para(:,y)), 100);
    y_fit = pdf(pd,x);
    histogram(para(:,y),100,'Normalization','pdf','FaceColor','b','EdgeColor','none');
    plot(x, y_fit, 'r', 'LineWidth', 2); % Plot the curve in red
    hold on
    xlabel('\gamma_t_1')
    ylabel('Frequency')
    z=z+1;
    y=y+1;
    subplot(n,2,z)
    plot(para(:,y));
    axis tight;
    xlabel('Iteration')
    ylabel('\gamma_t_2')
    hold on
    z=z+1;
    subplot(n,2,z)
    hold on
    pd = fitdist(para(:,y), 'Normal'); 
    x = linspace(min(para(:,y)), max(para(:,y)), 100);
    y_fit = pdf(pd,x);
    histogram(para(:,y),100,'Normalization','pdf','FaceColor','b','EdgeColor','none');
    plot(x, y_fit, 'r', 'LineWidth', 2); % Plot the curve in red
    hold on
    % xlabel('\gamma_t_2')
    ylabel('Frequency')
    z=z+1;
    y=y+1;
end

if(strcmp(Model_type(2,1),'Stationary'))
   
    subplot(n,2,z)
    plot(para(:,y));
    axis tight;
    xlabel('Iteration')
    ylabel('\sigma')
    hold on
    z=z+1;
    subplot(n,2,z)
    hold on
    pd = fitdist(para(:,y), 'Normal'); 
    x = linspace(min(para(:,y)), max(para(:,y)), 100);
    y_fit = pdf(pd,x);
    histogram(para(:,2),100,'Normalization','pdf','FaceColor','b','EdgeColor','none');
    plot(x, y_fit, 'r', 'LineWidth', 2); % Plot the curve in red
    hold on
    xlabel('\sigma')
    ylabel('Frequency')
    z=z+1;
    y=y+1;
end

if(strcmp(Model_type(2,1),'Non-Stationary Trend'))
    param_names = cell(1, size(int_var,1));
for param_idx = 1:size(int_var,1)
    param_names{param_idx} = ['\sigma_t_', num2str(param_idx)];
end
    for i=1:1:size(int_var,1)
    
    
    subplot(n,2,z)
    plot(para(:,y));
    axis tight;
    xlabel('Iteration')
    ylabel(param_names{i})
    hold on
    z=z+1;
    subplot(n,2,z)
    hold on
    pd = fitdist(para(:,y), 'Normal'); 
    x = linspace(min(para(:,y)), max(para(:,y)), 100);
    y_fit = pdf(pd,x);
    histogram(para(:,y),100,'Normalization','pdf','FaceColor','b','EdgeColor','none');
    plot(x, y_fit, 'r', 'LineWidth', 2); % Plot the curve in red
    hold on
    % xlabel(param_names{i})
    ylabel('Frequency')
    z=z+1;
    y=y+1;
    end
end

if(strcmp(Model_type(3,1),'Stationary'))
    subplot(n,2,z)
    plot(para(:,y));
    axis tight;
    xlabel('Iteration')
    ylabel('\mu')
    hold on
    z=z+1;
    subplot(n,2,z)
    hold on
    pd = fitdist(para(:,y), 'Normal'); 
    x = linspace(min(para(:,y)), max(para(:,y)), 100);
    y_fit = pdf(pd,x);
    histogram(para(:,2),100,'Normalization','pdf','FaceColor','b','EdgeColor','none');
    plot(x, y_fit, 'r', 'LineWidth', 2); % Plot the curve in red
    hold on
    % xlabel('\mu')
    ylabel('Frequency')
    z=z+1;
    y=y+1;
end

if(strcmp(Model_type(3,1),'Non-Stationary Trend'))
    param_names = cell(1, size(int_var,1));
for param_idx = 1:size(int_var,1)
    param_names{param_idx} = ['\mu_t_', num2str(param_idx)];
end
for i=1:1:size(int_var,1)
    subplot(n,2,z)
    plot(para(:,y));
    axis tight;
    xlabel('Iteration')
    ylabel(param_names{i})
    hold on
    z=z+1;
    subplot(n,2,z)
    hold on
    pd = fitdist(para(:,y), 'Normal'); 
    x = linspace(min(para(:,y)), max(para(:,y)), 100);
    y_fit = pdf(pd,x);
    histogram(para(:,y),100,'Normalization','pdf','FaceColor','b','EdgeColor','none');
    plot(x, y_fit, 'r', 'LineWidth', 2); % Plot the curve in red
    hold on
    % xlabel(param_names{i})
    ylabel('Frequency')
    z=z+1;
    y=y+1;

end
end
end
case 'gp'
    