% An example for the Bayesian variable selection method in AFT model for survival data,  
% developed in Zhang, Z., Sinha, S., Maiti, T., and Shipp, E. (2016). Bayesian variable selection  
% in the AFT model with an application to the SEER breast cancer data. To appear in the  
% Statistical Methods in Medical Research.  The proposed method handles the case where the number of 
% parameters is fixed and is smaller than the sample size n.     
% Please contact the authors if there are any questions or implementation issues:  
%    Zhen Zhang, zhangz19@galton.uchicago.edu
%
%
%
%
%
% In this file we illustrate how to use the main Matlab function AFT_Bayes_LASSO.
%
%
%======== The input variables of this function are 
% V: observed time to event, in terms of notation, V=min(T, C), where T and C are time-to-event and the censoring
%                   time, respectively.
% Delta: binary vector indicating if an observation is not censored, i.e. Delta = (T<=C);
% Z: n ny q design matrix with q candidate variables
% N: truncated maximal number of clusters for Dirichlet Process prior
% tot: number of total iterations per MCMC run
% burn: number of iterations as burn-in time
% init_beta: q by 1 initial value of the Beta-parameter, can be obtained
%        from parametric AFT (see R package flexsurv)
% randomSeed: a single number for the seed of random number generator
%
%
%
%======== The output variables of this function are 
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

id = 1; % Indexing different replicates of simulated data
rng('default'); rng(8*id)

% Simulate data with scaled log-Weibull distributed error. 
% Similar to scenario 1 in Zhang et al. (2016). See the article for details
n = 5e3;
q = 10;
p = 4; 
inds = 1:p;
betas = zeros(q,1); betas(inds) = [0.5, 0.5, 0.35, -0.35]';
Z = binornd(1, 0.5, [n,q]);
mu = Z*betas;
labs = ones(1,n);
lambda0 = 29.2672; k0 = 0.8178; C0 = 1.57;
eps = log(wblrnd(lambda0, k0,[n,1]))/C0;
eps = exp(-0.5*mu.^2).*eps;
K = 78;
T = exp(1 + mu + eps);
C = Z(:,1) + Z(:,2) + unifrnd(0,K,[n,1]);
V = min(T,C);
Delta = (T<=C);
fprintf('The proportion of non-censoring: %3.4f\n', mean(Delta))

% check the empirical survival rates for some example groups
n0 = 4;
Z0 = zeros(n0, q); 
Z0(1, 1:p) = [1,1,1,0]; % survival group with the highest survival rate
Z0(2, 1:p) = [1,0,1,1]; % survival group with an intermediate survival rate
Z0(3, 1:p) = [1,1,1,1]; % survival group with an intermediate survival rate
Z0(4, 1:p) = [0,0,0,1]; % survival group with the worst survival rate

t0 = 0:0.2:60; len = numel(t0); trueY = zeros(len, n0);
for i = 1:n0
    trueMu = Z0(i,:)*betas;
    x = exp(0.5*trueMu^2)*(log(t0)-1-trueMu); 
    trueY(:,i) = exp(- (exp(C0*x)/lambda0).^k0);
end
Trues.eps = eps; Trues.T = T; Trues.C = C;
Trues.beta = betas; Trues.inds = inds;
Trues.labs = labs;

% store results
tot = 200; burn = 100; % for pilot run
% %========== in practice, need to modify total number of iterations per chain, and burn-in period
% tot = 2e4; burn = 15e3; 

N =150; %modify the truncated maximal number of clusters for Dirichlet Process prior 
nch = 3; % number of MCMC chains
myseeds = [1,2,3]; %set seeds for random number generator for each chain
ngammas = 2; 
niter = tot-burn; 
nsample = niter*nch; 
matpara = zeros(nsample, q*2+1+ngammas+1 + 1);
Ps = zeros(nsample, N);
Theta = zeros(nsample, N*2);
initBetas = Trues.beta; 

% run MCMC chains, may be time-consuming with supplied "tot" values.
for ch = 1:nch
    i0 = (ch-1)*niter + (1:niter); 
    fprintf('\nMCMC chain %d: ',ch)
    [matpara(i0,:), Theta(i0,:), Ps(i0,:)] = AFT_Bayes_LASSO(V, Delta, Z, N, tot, burn, initBetas, myseeds(ch)); 
end

% check beta-parameter coverage
betaSamples = matpara(:, 1:q)'; 
lb = quantile(betaSamples,0.025,2); ub = quantile(betaSamples,0.975,2); 
fprintf('\nTrue, posterior mean, lower and upper bound of 95%% credible interval, coverage of the true:\n')
disp(num2str([Trues.beta,mean(betaSamples,2), lb, ub, (lb<Trues.beta).*(ub>Trues.beta)==1], 3))

% check survival probability
mu0 = Z0*betaSamples;
nu0 = zeros(n0, nsample);
spm = nan(len, n0); spl = nan(len, n0); spu = nan(len, n0);
fprintf('\nThe design matrix for some example groups for survival probability:\n')
disp(Z0)

for i = 1:n0
    for j = 1:nsample
        nu0(i, j) = exp( [mu0(i,j), mu0(i,j).^2]*matpara(j, (2*q+1)+(1:ngammas))' );
    end
    for t = 1:len
        t1 = t0(t);
        tmp = 1 - sum(Ps.*normcdf(...
            ( repmat( ((log(t1) - mu0(i,:))./sqrt(nu0(i,:)))', [1,N]) - Theta(:,1:N) )./sqrt(Theta(:,N+(1:N)))...
            ), 2);
        spm(t, i) = mean(tmp);  spl(t, i) = quantile(tmp, .025); spu(t, i) = quantile(tmp, .975);
    end
    
    % Plot the true, empirical, and estimated survival probability
    subplot(2,2,i); plot(t0, trueY(:,i), 'r-')
    xlim([Z0(i,:)*betas+1, 60])
    hold on,
    for j = t0
        plot(j, mean(V(~( sum((Z(:,1:p) - repmat(Z0(i,1:p), [n,1])).^2, 2) ))>=j), 'k.')
    end
    plot(t0, spm(:,i), 'b-')
    plot(t0, spl(:,i), 'b--')
    plot(t0, spu(:,i), 'b--')
    hold off
end

% % save the figure
% fac = .8;
% set(gcf, 'PaperPosition', [0 0 4 4]/fac);
% set(gcf, 'PaperSize', [4 4]/fac);
% set(gca,'FontSize',12)
% saveas(gcf, strcat('./','fig',num2str(id)), 'pdf')

disp('Demo completed. For full results, need to increase tot for achieving reasonable mixing rates and convergence')
