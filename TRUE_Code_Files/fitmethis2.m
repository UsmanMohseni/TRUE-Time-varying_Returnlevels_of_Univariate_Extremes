function F = fitmethis2(data, varargin)
% FITMETHIS finds best-fitting distribution and includes KS test results.

% Defaults & Storage
dtype = 'cont';
ntrials = [];
fig = 'off';
alpha = 0.05;
criterion = 'KS';
output = 'off';
F = struct('name', {}, 'par', [], 'ci', [], 'LL', [], 'aic', [], 'ksstat', [], 'kspval', []);
prefdist = [];
plotdata = 'off';
plotdist = 1;
applykernel = false;

% Arguments
for j = 1:2:length(varargin)
    string = lower(varargin{j});
    switch string(1:min(4,length(string)))
        case 'dtyp'
            dtype = varargin{j+1};
        case 'ntri'
            ntrials = varargin{j+1};
        case 'figu'
            fig = varargin{j+1};
        case 'alph'
            alpha = varargin{j+1};
        case 'crit'
            criterion = varargin{j+1};
        case 'outp'
            output = varargin{j+1};
        case 'pref'
            prefdist = varargin{j+1};
        case 'pdat'
            plotdata = varargin{j+1};
        case 'pdis'
            plotdist = varargin{j+1};
        case 'kern'
            applykernel = varargin{j+1};
        otherwise
            error('Unknown argument name');
    end
end

% Distributions
Cdist = {'gamma'; 'gev';'gp'; 'lognormal'};
mustbepos = 11;

% Fit all... 
switch dtype(1:4)

    % Continuous
    case 'cont'
        for j = 1:numel(Cdist)

            % If negative values, only fit normal
            if min(data) < 0
                [phat, pci] = mle(data, 'distribution', 'normal', 'alpha', alpha);
                F(j).name = Cdist{j};
                F(j).par = phat;
                F(j).ci = pci;
                pdfv = pdf('normal', data, F(j).par(1), F(j).par(2));
                F(j).LL = sum(log(pdfv(pdfv > 0 & ~isinf(pdfv))));
                F(j).aic = 2*2 - 2*F(j).LL;
                [F(j).ksstat, F(j).kspval] = kstest((data - mean(data)) / std(data)); % KS test
                break

            % Check: if values > 1 for Beta, do nothing
            elseif strcmp('beta', Cdist{j}) && max(data) > 1
                F(j).name = 'beta';
                F(j).LL = -Inf;
                F(j).aic = Inf;
                F(j).ksstat = NaN;
                F(j).kspval = NaN;

            % Check: if values > 0 for some distr. (they are sorted), do nothing
            elseif j >= mustbepos && min(data) == 0
                F(j).name = Cdist{j};
                F(j).LL = -Inf;
                F(j).aic = Inf;
                F(j).ksstat = NaN;
                F(j).kspval = NaN;

            % Any other case do the fit ...
            else
                try
                    [phat, pci] = mle(data, 'distribution', Cdist{j}, 'alpha', alpha);
                    F(j).name = Cdist{j};
                    F(j).par = phat;
                    F(j).ci = pci;
                    if numel(F(j).par) == 1
                        pdfv = pdf(F(j).name, data, F(j).par(1));
                    elseif numel(F(j).par) == 2
                        pdfv = pdf(F(j).name, data, F(j).par(1), F(j).par(2));
                    else
                        pdfv = pdf(F(j).name, data, F(j).par(1), F(j).par(2), F(j).par(3));
                    end
                    F(j).LL = sum(log(pdfv(pdfv > 0 & ~isinf(pdfv))));
                    F(j).aic = 2 * numel(F(j).par) - 2 * F(j).LL;

                    % KS test
                    % [F(j).ksstat, F(j).kspval] = kstest((data - mean(data)) / std(data));
                switch Cdist{j}
                    case 'gamma'
                        params=gamfit(data);
                        cdf=gamcdf(data,params(1),params(2));
                    [~,F(j).kspval, F(j).ksstat,~] = kstest(data,"CDF",[data cdf]);
                    case 'gev'
                        params=gevfit(data);
                        cdf=gevcdf(data,params(1),params(2),params(3));
                    [~,F(j).kspval, F(j).ksstat,~] = kstest(data,"CDF",[data cdf]);
                    case 'gp'
                        params=gpfit(data);
                        cdf=gpcdf(data,params(1),params(2));
                    [~,F(j).kspval, F(j).ksstat,~] = kstest(data,"CDF",[data cdf]);
                    case 'lognormal'
                            params=lognfit(data);
                        cdf=logncdf(data,params(1),params(2));
                    [~,F(j).kspval, F(j).ksstat,~] = kstest(data,"CDF",[data cdf]);
                end
                catch
                    F(j).name = Cdist{j};
                    F(j).par = NaN;
                    F(j).LL = -Inf;
                    F(j).aic = Inf;
                    F(j).ksstat = NaN;
                    F(j).kspval = NaN;
                    fprintf('%s%s\n', Cdist{j}, ' distribution not applicable')
                end
            end

        end

end

% Order by criterion
switch criterion
    case 'LL'
        index = sortrows([(1:size(F,2))', [F.LL]'], -2);
    case 'AIC'
        index = sortrows([(1:size(F,2))', [F.aic]'], 2);
    case 'KS'
        index = sortrows([(1:size(F,2))', [F.ksstat]'], 2);
end
F = F(index(:,1));

% Nice screen output
if strcmp('on', output)
    fprintf('\n\t\t\t\tName\t\tPar1\t\tPar2\t\tPar3\t\tLogL\t\tAIC\t\tKS Stat\t\tKS p-value\n')
    for j = 1:size(F,2)
        switch numel(F(j).par)
            case 1
                fprintf('%20s \t%10.3e \t\t\t\t\t\t\t%10.3e \t%10.3e \t%10.3e \t%10.3e\n', ...
                    F(j).name, F(j).par, F(j).LL, F(j).aic, F(j).ksstat, F(j).kspval)
            case 2
                fprintf('%20s \t%10.3e \t%10.3e \t\t\t\t%10.3e \t%10.3e \t%10.3e \t%10.3e\n', ...
                    F(j).name, F(j).par, F(j).LL, F(j).aic, F(j).ksstat, F(j).kspval)
            case 3
                fprintf('%20s \t%10.3e \t%10.3e \t%10.3e \t%10.3e \t%10.3e \t%10.3e \t%10.3e\n', ...
                    F(j).name, F(j).par, F(j).LL, F(j).aic, F(j).ksstat, F(j).kspval)
        end
    end
end

% End of fitmethis
end
