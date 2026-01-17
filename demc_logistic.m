%%  Function for DE-MC (LOGISTIC)
function [para,Rhat,accrate,mix] = demc_logistic( ...
    Model_type,int_var,AMS,evl,bur,sts,cha,Covariate,int_cp_index)

% ==========================================================
% Determine number of parameters (dim)
% ==========================================================
[m,n] = size(int_var);
dim = 0;
for i = 1:m
    for j = 1:n
        if ~isnan(int_var(i,j))
            dim = dim + 1;
        end
    end
end

% ==========================================================
% Priors: user-defined or defaults
% ==========================================================
if evalin('base',"exist('Logistic_Settings_Matrix','var')") && ...
        ~isempty(evalin('base','Logistic_Settings_Matrix'))

    Logistic_Settings_Matrix = evalin('base','Logistic_Settings_Matrix');

    % Location
    loc_dist = Logistic_Settings_Matrix{1,1};
    loc_p1   = Logistic_Settings_Matrix{2,1};
    loc_p2   = Logistic_Settings_Matrix{3,1};

    % Scale
    scale_dist = Logistic_Settings_Matrix{1,2};
    scale_p1   = Logistic_Settings_Matrix{2,2};
    scale_p2   = Logistic_Settings_Matrix{3,2};
else
    loc_dist   = 'norm';
    loc_p1     = 0;
    loc_p2     = 1000;
    scale_dist = 'norm';
    scale_p1   = 0;
    scale_p2   = 1000;
end

% ==========================================================
% Logistic log-PDF (numerically stable)
% ==========================================================
logistic_logpdf = @(x,mu,s) ...
    -log(s) - (x-mu)./s - 2*log1p(exp(-(x-mu)./s));

% ==========================================================
% DE-MC initialization
% ==========================================================
accept = zeros(evl,1);
reject = zeros(evl,1);
mix    = zeros(cha,dim+1,evl);
z_z    = zeros(cha,dim);

obs_int = repmat(AMS',cha,1);

% ==========================================================
% Initial chains (priors)
% ==========================================================
[p_x_k,mix_k,x_k] = Prior_estimator( ...
    int_var(:,1:2),int_cp_index,AMS,Covariate,obs_int, ...
    1000,cha,Model_type{1,1},loc_dist,loc_p1,loc_p2);

[p_x_sigma,mix_x_sig,x_sigma] = Prior_estimator( ...
    int_var(:,3:4),int_cp_index,AMS,Covariate,obs_int, ...
    1000,cha,Model_type{2,1},scale_dist,scale_p1,scale_p2);

% ==========================================================
% Initial likelihood (LOGISTIC)
% ==========================================================
non_xx = logistic_logpdf(obs_int,x_k,exp(x_sigma));
likeli_joint = sum(non_xx,2);

p_post = likeli_joint + p_x_sigma + p_x_k;

mix(:,1:dim,1) = [mix_k mix_x_sig];
mix(:,end,1)   = p_post;

JR(:,1) = 2.38 ./ sqrt(2*(1:dim)');

% ==========================================================
% DE-MC Sampling
% ==========================================================
for i = 1:evl

    [~,ttt] = sort(rand(cha-1,cha));
    D = rand(cha,dim);

    for jj = 1:cha
        ii = ones(cha,1); ii(jj)=0;
        idx = find(ii>0);
        rr  = idx(ttt(1:2,jj));

        c = find(D(jj,1:dim)>0);
        if isempty(c)
            c = randperm(dim,1);
        end
        Dim = numel(c);

        if rand < 0.8
            g = JR(Dim);
        else
            g = 1;
        end

        delta(jj,c) = g*(mix(rr(1),1:dim,i) - mix(rr(2),1:dim,i));

        if sum(delta(jj,:).^2) == 0
            [RR,P] = chol(cov(mix(:,1:dim,i)) + 1e-5*eye(dim));
            if P == 0
                delta(jj,:) = randn(1,dim)*(2.38/sqrt(dim))*RR;
            end
        end
    end

    z_z(:,1:dim) = mix(:,1:dim,i) + delta(:,1:dim);
    p_post = mix(:,end,i);

    % ======================================================
    % Prior simulation
    % ======================================================
    if strcmp(Model_type{1,1},'Stationary')
        [p_z_k,z_k] = Prior_Simulation( ...
            int_var(:,1:2),int_cp_index,AMS,Covariate,obs_int, ...
            1000,cha,z_z(:,1),Model_type{1,1},loc_dist,loc_p1,loc_p2);
    else
        [p_z_k,z_k] = Prior_Simulation( ...
            int_var(:,1:2),int_cp_index,AMS,Covariate,obs_int, ...
            1000,cha,z_z(:,1:m),Model_type{1,1},loc_dist,loc_p1,loc_p2);
    end

    [p_z_sigma,z_sigma] = Prior_Simulation( ...
        int_var(:,3:4),int_cp_index,AMS,Covariate,obs_int, ...
        1000,cha,z_z(:,m+1:end),Model_type{2,1}, ...
        scale_dist,scale_p1,scale_p2);

    % ======================================================
    % Logistic likelihood
    % ======================================================
    non_zxx = logistic_logpdf(obs_int,z_k,exp(z_sigma));
    p_likeli_joint = sum(non_zxx,2);

    pz_post = p_likeli_joint + p_z_sigma + p_z_k;

    ratio = pz_post - p_post;
    alfa  = min(exp(ratio),1);

    r = rand(cha,1);
    cidx = find(alfa > r);
    ridx = find(alfa <= r);

    mix(cidx,:,i+1) = [z_z(cidx,:), pz_post(cidx)];
    mix(ridx,:,i+1) = mix(ridx,:,i);

    accept(i) = numel(cidx);
    reject(i) = numel(ridx);
end

% ==========================================================
% Posterior extraction
% ==========================================================
for i = 1:dim
    P = mix(:,i,bur:evl);
    para(:,i) = reshape(P,1,[]).';
    para(:,i) = para(sts:end,i);
end

sur = 'nonsta';
Rhat = Rconver(1,evl,bur,cha,mix,sts,sur,dim);
accrate = sum(accept) ./ sum(accept+reject);
disp(accrate)

end
