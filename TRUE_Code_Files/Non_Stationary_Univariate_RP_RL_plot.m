function Non_Stationary_Univariate_RP_RL_plot(params,best_fit,Time,Label,no)

    X = ceil(prctile(Time, 95)); % taking 95th percentile value as design parameters
     RGB = [0.93, 0.86, 0.55;  % Soft yellow
           0.2, 0.4, 0.6;    % Muted blue
           0.8, 0.4, 0.4];   % Muted red

           
    switch best_fit
        case 'gev'
            shape = params{1};
            scale_t = params{2};
            mu_t = params{3};
            scale_t = scale_t(:, X);
            shape = shape(:, X);
            mu_t = mu_t(:, X);
            para = [shape exp(scale_t) mu_t];
            z = norminv(1 - 0.05 / 2);
            CI(1:3, 1:2) = NaN;
            for i = 1:1:3
                mu(i, 1) = mean(para(:, i));
                sd = std(para(:, i));
                CI(i, 1) = mu(i, 1) - z * sd;
                CI(i, 2) = mu(i, 1) + z * sd;
            end
            z2 = norminv(1 - 0.01 / 2);
            CI2(1:3, 1:2) = NaN;
            for i = 1:1:3
                mu(i, 1) = mean(para(:, i));
                sd = std(para(:, i));
                CI2(i, 1) = mu(i, 1) - z2 * sd;
                CI2(i, 2) = mu(i, 1) + z2 * sd;
            end
            for i = 2:1:100
                RL1(i - 1, :) = (gevinv(1 - 1 / i, mu(1, 1), mu(2, 1), mu(3, 1)));
                RL2(i - 1, :) = (gevinv(1 - 1 / i, CI(1, 1), CI(2, 1), CI(3, 1)));
                RL3(i - 1, :) = (gevinv(1 - 1 / i, CI(1, 2), CI(2, 2), CI(3, 2)));
            end
            for i = 2:1:100
                RL4(i - 1, :) = (gevinv(1 - 1 / i, CI2(1, 1), CI2(2, 1), CI2(3, 1)));
                RL5(i - 1, :) = (gevinv(1 - 1 / i, CI2(1, 2), CI2(2, 2), CI2(3, 2)));
            end
        case 'gp'
            shape = params{1};
            scale_t = params{2};
            scale_t = scale_t(:, X);
            shape = shape(:, X);
            para = [exp(shape) scale_t];
            z = norminv(1 - 0.05 / 2);
            CI(1:2, 1:2) = NaN;
            for i = 1:1:2
                mu(i, 1) = mean(para(:, i));
                sd = std(para(:, i));
                CI(i, 1) = mu(i, 1) - z * sd;
                CI(i, 2) = mu(i, 1) + z * sd;
            end
            z2 = norminv(1 - 0.01 / 2);
            CI2(1:3, 1:2) = NaN;
            for i = 1:1:2
                mu(i, 1) = mean(para(:, i));
                sd = std(para(:, i));
                CI2(i, 1) = mu(i, 1) - z2 * sd;
                CI2(i, 2) = mu(i, 1) + z2 * sd;
            end
            for i = 2:1:100
                RL1(i - 1, :) = (gpinv(1 - 1 / i, mu(1, 1), mu(2, 1)));
                RL2(i - 1, :) = (gpinv(1 - 1 / i, CI(1, 1), CI(2, 1)));
                RL3(i - 1, :) = (gpinv(1 - 1 / i, CI(1, 2), CI(2, 2)));
            end
            for i = 2:1:100
                RL4(i - 1, :) = (gpinv(1 - 1 / i, CI2(1, 1), CI2(2, 1)));
                RL5(i - 1, :) = (gpinv(1 - 1 / i, CI2(1, 2), CI2(2, 2)));
            end
        case 'lognormal'
            shape = params{1};
            scale_t = params{2};
            scale_t = scale_t(:, X);
            shape = shape(:, X);
            para = [(shape) exp(scale_t)];
            z = norminv(1 - 0.05 / 2);
            CI(1:2, 1:2) = NaN;
            for i = 1:1:2
                mu(i, 1) = mean(para(:, i));
                sd = std(para(:, i));
                CI(i, 1) = mu(i, 1) - z * sd;
                CI(i, 2) = mu(i, 1) + z * sd;
            end
            z2 = norminv(1 - 0.01 / 2);
            CI2(1:3, 1:2) = NaN;
            for i = 1:1:2
                mu(i, 1) = mean(para(:, i));
                sd = std(para(:, i));
                CI2(i, 1) = mu(i, 1) - z2 * sd;
                CI2(i, 2) = mu(i, 1) + z2 * sd;
            end
            for i = 2:1:100
                RL1(i - 1, :) = (logninv(1 - 1 / i, mu(1, 1), mu(2, 1)));
                RL2(i - 1, :) = (logninv(1 - 1 / i, CI(1, 1), CI(2, 1)));
                RL3(i - 1, :) = (logninv(1 - 1 / i, CI(1, 2), CI(2, 2)));
            end
            for i = 2:1:100
                RL4(i - 1, :) = (logninv(1 - 1 / i, CI2(1, 1), CI2(2, 1)));
                RL5(i - 1, :) = (logninv(1 - 1 / i, CI2(1, 2), CI2(2, 2)));
            end
        case 'gamma'
            shape = params{1};
            scale_t = params{2};
            scale_t = scale_t(:, X);
            shape = shape(:, X);
            para = [exp(shape) exp(scale_t)];
            z = norminv(1 - 0.05 / 2);
            CI(1:2, 1:2) = NaN;
            for i = 1:1:2
                mu(i, 1) = mean(para(:, i));
                sd = std(para(:, i));
                CI(i, 1) = mu(i, 1) - z * sd;
                CI(i, 2) = mu(i, 1) + z * sd;
            end
            z2 = norminv(1 - 0.01 / 2);
            CI2(1:3, 1:2) = NaN;
            for i = 1:1:2
                mu(i, 1) = mean(para(:, i));
                sd = std(para(:, i));
                CI2(i, 1) = mu(i, 1) - z2 * sd;
                CI2(i, 2) = mu(i, 1) + z2 * sd;
            end
            for i = 2:1:100
                RL1(i - 1, :) = (gaminv(1 - 1 / i, mu(1, 1), mu(2, 1)));
                RL2(i - 1, :) = (gaminv(1 - 1 / i, CI(1, 1), CI(2, 1)));
                RL3(i - 1, :) = (gaminv(1 - 1 / i, CI(1, 2), CI(2, 2)));
            end
            for i = 2:1:100
                RL4(i - 1, :) = (gaminv(1 - 1 / i, CI2(1, 1), CI2(2, 1)));
                RL5(i - 1, :) = (gaminv(1 - 1 / i, CI2(1, 2), CI2(2, 2)));
            end
    end

    % Plotting
    figure(no)

  
    X = (2:100)';
    % Fill for the confidence interval
    fill([X; flipud(X)], [RL4; flipud(RL5)], RGB(1,:), 'EdgeColor', 'none') % Soft yellow fill
    hold on
    % Median return level
    plot(X, RL1, 'Color', RGB(2,:), 'LineWidth', 3, 'DisplayName', 'Median') % Muted blue curve
    hold on
    % Confidence intervals
    plot(X, RL2, '-.', 'Color', RGB(3,:), 'LineWidth', 2, 'DisplayName', '95% CI Lower') % Muted red dashed
    hold on
    plot(X, RL3, '-.', 'Color', RGB(3,:), 'LineWidth', 2, 'DisplayName', '95% CI Upper') % Muted red dashed

    % Setting titles and labels
    title('Non-Stationary: Return Period vs Return Level', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
    xlabel('Return Period (years)', 'FontSize', 12, 'FontName', 'Times New Roman')
    ylabel(Label, 'FontSize', 12, 'FontName', 'Times New Roman')
    legend('Ensembles', 'Mean ', '95% CI', '5% CI')
    % Adding legend
    legend('Location', 'northwest', 'FontSize', 10, 'FontName', 'Times New Roman', 'Box', 'off')
    grid on
set(gcf, 'Color', 'w', 'Units', 'centimeters', 'Position', [18.12, 4.25, 17.75,7.6]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman') % Adjusting font for axes
end

