function [matpara, Theta, Ps] = AFT_Bayes_LASSO(V, Delta, Z, N, tot, burn, init_beta, randomSeed)
% MATLAB code for the Bayesian variable selection method in AFT model for survival data,
% developed in  Zhang, Z., Sinha, S., Maiti, T., and Shipp, E. (2016) Bayesian variable selection
% in the AFT model with an application to the SEER breast cancer data. To appear in the
% Statistical Methods in Medical Research.
% Please contact the authors if there are any questions or implementation issues:
%    Zhen Zhang, zhangz19@galton.uchicago.edu
% In addition to supplying the input variables, the hyperparameters may
% also need to be modified for incorporating the prior information
% ================== List of input variables:
%    --  V: n by 1 vector of the observed survival time with sample size n
%    --  Delta: n by 1 vector indicating if an observation is not censored
%    --  Z: n ny q design matrix with q candidate variables
%    --  N: truncated maximal number of clusters for Dirichlet Process prior
%    --  tot: number of total iterations per MCMC run
%    --  burn: number of iterations as burn-in time
%    --  init_beta: q by 1 initial value of the Beta-parameter, can be obtained
%                   from parametric AFT (see R package flexsurv)
%    --  randomSeed: a single number for the seed of random number generator
% ================== List of output variables:
% matpara: a matrix with K ( = tot-burn) rows that records MCMC samples in the order of:
%    -- q columns: Beta-parameter for the slopes of the q candidate variables
%    -- q columns: u-parameter for the precision of the q candidate variables
%    -- 1 column: lambda-parameter for the LASSO penalty
%    -- ngammas columns: gamma-parameter for the variance modeling.
%                  Default ngammas = 2 for quadratic function
%    -- 2 columns: alpha and zeta parameter in the paper
% Theta: K by 2*N matrix that records mean \theta_1 (first N
%          columns) and variance \theta_2 (last N columns) for the N clusters
% Ps: K by N matrix that records probability of each DP cluster

% set user-specified hyperparameters
invtau2_eta = 1/1e-3; % precision for eta prior, i.e. 1/\tau_{\eta}^2
phi1 = 0; phi2 = 2.5; phi3 = 1; % specify \phi_1, \phi_2 and \phi_3 in the paper
mean_IG = 1; var_IG = .01; % specify prior mean and variance of \zeta (Inverse-Gamma prior)
phi4 = 2+mean_IG^2/var_IG;
phi5 = mean_IG*(phi4-1);
mean_gamma = 0; var_gamma2 = 1e4; % specify prior mean and variance of \gamma
lambda_mean = 0.1; lambda_var = 5; % specify prior mean and variance of \lambda^2 (Gamma prior)
b_lambda = lambda_var/lambda_mean;
a_lambda = lambda_mean/b_lambda;
b_lambda = 1/b_lambda;
alpha_mean = 1; alpha_var = 10; % specify prior mean and variance of \alpha (Gamma prior)
b_alpha = alpha_var/alpha_mean;
a_alpha = alpha_mean/b_alpha;

verbose = 0; %track MCMC samples
reportrate = 0;
ngammas = 2;
[n,q] = size(Z);
cind = find(Delta==0);
L_cind = length(cind);

rng('default'); rng(randomSeed)

% initialize parameters
lTstar = log( V + unifrnd(0,1,[n,1]).*(Delta==0) );
alphas = 10;
lambda2 = gamrnd(a_lambda, 1/b_lambda);
eps_eta = normrnd(0,sqrt(1/invtau2_eta),[n,1]);
betas = init_beta; %(Z'*Z)\(Z'*lTstar); % may replace with more reasonable guess, say, regular AFT estimates
gammas = zeros(1,ngammas);
ps = 1/N*ones(1,N);
labs = ones(1,n); %randsample(N, n, true, ps);
u = ones(1,q);
zeta = 1;
theta = zeros(N,2);
theta(:,2) = 1./gamrnd(phi2, 1/phi3, [N,1]);
theta(:,1) = normrnd(phi1*ones(N,1), sqrt(zeta*theta(:,2)));
mu = Z*betas;
eta = mu + eps_eta;
nu = exp([eta,eta.^2]*gammas');
fulltheta = theta(labs, 1);
fulltheta2 = 1./( nu.*theta(labs, 2) );

% Adaptive MCMC set-up
objrate = 0.44;
batchLen = min(50, tot); batchNum = 0;
batchTot = tot/batchLen;
nlen = ngammas + n; %gammas, n eta
accepts = zeros(1, nlen);
rates = zeros(batchTot, nlen);
% set initial log- standard deviation of the proposal density
tunings = [-4*ones(1,ngammas), -0.5*ones(1,n)]; % may need to modify for reasonable initial acceptance rates

% store results
matpara = zeros(tot-burn, q*2+1+ngammas+1 + 1);
Ps = zeros(tot-burn, N);
Theta = zeros(tot-burn, N*2);

% MCMC run
tic
for iter = 1:tot
    
    if verbose == 1
        fprintf('%6d', iter)
        if(~mod(iter,20))
            fprintf('\n')
        end
    end
    
    % step 1 & 2: sample beta and u
    Ztild = repmat(sqrt(fulltheta2+invtau2_eta),[1,q]).*Z;
    Sigma = diag(u) + Ztild'*Ztild;
    Sigma = chol(Sigma, 'lower');
    Mu =  Sigma\( Z'*((lTstar - fulltheta.*sqrt(nu)).*fulltheta2 + eta*invtau2_eta) );
    betas = Sigma'\( randn(size(Mu)) + Mu );
    betas = betas';
    for j = 1:q
        u0 = rand(1); u_mu = sqrt(lambda2)/abs(betas(j));
        y0 = (randn(1))^2; tlambda2 = 2*lambda2;
        unew = u_mu + u_mu^2*y0/tlambda2 - u_mu/tlambda2*sqrt(2*tlambda2*u_mu*y0 + u_mu^2*y0^2);
        if u0 <= u_mu/(u_mu+unew)
            u(j) = unew;
        else
            u(j) = u_mu^2/unew;
        end
    end
    
    % step3: sample lambda2
    lambda2 = gamrnd(a_lambda+q, 1./(b_lambda + 0.5*sum(u.^2)));
    
    % step 4: sample gamma
    err = lTstar - Z*betas';
    fulltheta2 = 1./theta(labs, 2);
    for j = 1:ngammas
        gamma1 = gammas; gamma1(j) = normrnd(gammas(j), exp(tunings(j)));
        lognu = [eta,eta.^2]*gammas';
        lognu1 = [eta,eta.^2]*gamma1';
        loglik = - 0.5*sum((err./sqrt(exp(lognu)) - fulltheta).^2.*fulltheta2) - 0.5*sum(lognu) - (gammas(j)-mean_gamma)^2./var_gamma2;
        loglik1 = - 0.5*sum((err./sqrt(exp(lognu1)) - fulltheta).^2.*fulltheta2) - 0.5*sum(lognu1) - (gamma1(j)-mean_gamma)^2./var_gamma2;
        MH = exp(loglik1 - loglik);
        u0 = rand(1);
        if u0 <= MH
            gammas = gamma1; lognu = lognu1;
            accepts(j) = accepts(j)+1;
        end
    end
    nu = exp(lognu);
    
    % step 5: sample censored Ttar
    if ~isempty(cind)
        snu = sqrt(nu(cind));
        R = rand([L_cind,1]); tmp =  Z(cind,:)*betas' + snu.*fulltheta(cind);
        snu = snu./sqrt(fulltheta2(cind));
        lTstar(cind) = tmp + snu.*norminv( (1-R).*min(0.9999, normcdf((log(V(cind)) - tmp)./snu)) + R );
    end
    
    % step 6: sample eta
    mu = Z*betas';
    err = lTstar - mu;
    lognu = [eta, eta.^2]*gammas';
    eta1 = normrnd(eta, exp(tunings(ngammas+(1:n))'));
    lognu1 = [eta1, eta1.^2]*gammas';
    loglik = - 0.5*(err./sqrt(exp(lognu)) - fulltheta).^2.*fulltheta2 - lognu/2 - 0.5*invtau2_eta*(eta-mu).^2;
    loglik1 = - 0.5*(err./sqrt(exp(lognu1)) - fulltheta).^2.*fulltheta2 - lognu1/2 - 0.5*invtau2_eta*(eta1-mu).^2;
    MH = exp(loglik1 - loglik);
    u0 = rand(n,1); inds = find(u0 <= MH);
    if numel(inds) > 0
        eta(inds) = eta1(inds);
        lognu(inds) = lognu1(inds);
        accepts(ngammas+inds) = accepts(ngammas+inds)+1;
    end
    nu = exp(lognu);
    err = err./sqrt(nu);
    
    % step 7: sample labs
    ps0 = repmat(ps, [n,1]).*normpdf( (repmat(err, [1,N]) - repmat(theta(:,1)', [n,1]))./repmat(sqrt(theta(:,2)'), [n,1]) );
    ps0 = cumsum(ps0./repmat(sum(ps0,2), [1,N]), 2);
    myu0 = rand([n,1]);
    labs = sum(ps0 <= repmat(myu0, [1,N]), 2) + 1;
    ms = histc(labs, 1:N);
    
    % step 8: sample ps
    ps = gamrnd(alphas/N + ms, 1);
    ps = ps'./sum(ps);
    
    % step 9: sample thetas
    tmpsum = 0.0;
    for j = 1:N
        if ms(j) > 0
            var_theta = zeta/(1+ms(j)*zeta);
            mean_theta = var_theta*(phi1/zeta + sum(err(labs ==j)) );
            var_theta = theta(j,2)*var_theta;
            theta(j,1) = normrnd(mean_theta, sqrt(var_theta));
            theta(j,2) = 1./gamrnd(phi2 + 0.5*1   +  0.5*ms(j), ...
                1/(phi3 + 0.5*(theta(j,1)-phi1)^2/zeta + 0.5*sum((err(labs==j) - theta(j,1)).^2) ) );
            tmpsum = tmpsum + 0.5*(theta(j,1)-phi1)^2/theta(j,2);
        else % generate from base prob measure
            theta(j,2) = 1./gamrnd(phi2, 1/phi3 );
            theta(j,1) = normrnd(phi1, sqrt(zeta*theta(j,2)));
        end
    end
    zeta = 1/gamrnd(phi4 + 0.5*sum(ms>0), 1/(phi5 + tmpsum));
    fulltheta = theta(labs,1);
    fulltheta2 = 1./( nu.*theta(labs, 2) );
    
    % step 10: sample alpha
    alpha1 = gamrnd(a_alpha, b_alpha);
    sps = sum(log(ps));
    loglik =gammaln(alphas) - N*gammaln(alphas/N) +(alphas/N-1)*sps;
    loglik1 =gammaln(alpha1) - N*gammaln(alpha1/N) +(alpha1/N-1)*sps;
    MH = exp(loglik1 - loglik);
    u0 = rand(1);
    if u0 <= MH
        alphas = alpha1;
    end
    
    % record results after burn-in period
    if iter > burn
        iter0 = iter-burn;
        matpara(iter0,:) = [betas, u, lambda2, gammas, alphas, zeta];
        Theta(iter0,:) = reshape(theta, [1, 2*N]);
        Ps(iter0,:) = ps;
    end
    
    % tunning for adaptive MCMC
    if ~mod(iter, batchLen)
        batchNum = batchNum+1;
        accepts = accepts./batchLen;
        rates(batchNum,:) = accepts;
        if reportrate == 1
            disp(num2str( [ min(accepts(1:ngammas)), max(accepts(1:ngammas)), ...
                min(accepts(ngammas+(1:n))), max(accepts(ngammas+(1:n))) ],   2))
        end
        tunings = tunings + sign((accepts>objrate)-0.5).*min(0.01, 1/sqrt(batchNum));
        accepts = zeros(1,nlen);
    end
    
end

runtime = toc/60;
fprintf('%d iterations are done with elapsed time %.2f minutes.\n', tot, runtime)

% save(yourfilename,'matpara','Theta','Ps')

end


