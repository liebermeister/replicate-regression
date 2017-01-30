%function [x_max_posterior, neg_log_posterior, x_min, x_max, x_cov_posterior, y_cov_posterior] = maximise_bayesian_posterior_L1(data_mean, data_std, prior_mean, prior_std, R)
%
% Minimise -log posterior =   || (y(x)-ydata) / sigmay ||_L1 
%                           + || (x - xprior) / sigmax ||_L1
%
% where data vector y and parameter vector x are linearly related y(x) = R * x
%
% data_mean,  data_std : (n_data x 1) colunm vectors
% prior_mean, prior_std: (n_pars x 1) colunm vectors
%
% NOTE THAT THE SOLUTION NEED NOT BE UNIQUE!!
% Using linprog (the default) tends to give a symmetric solution; cplex may yield a non-symmetric solution
%
% NOTE that the posterior covariance matrices are not exact,  but
% computed via a simple heuristics:
% Start from the optimum point. 
% Allow for an (additive) increase of the L1 norm by 1/2
% Determine, for each parameter, the minimum and maximum parameter value that can still
%   be realised under this constraint
% Take 1/2 [xmax_i-xmin_i] as a "posterior standard deviation" and ignore correlations
% Turn these "standard deviations" into a covariance matrix, and derive from it the
% Covariance matrix of the data values y

function [x_max_posterior, neg_log_posterior, x_min, x_max, x_cov_posterior, y_cov_posterior] = maximise__bayesian_posterior_L1(data_mean, data_std, prior_mean, prior_std, R)

% bound for parameter values (only  make the lp problem formally bounded)
% if parameters are on log scale, this is huge; otherwise, this value may have to be increased
x_threshold = 1000; 

solver = 'linprog';
%solver = 'cplex';

if ~exist('cplexlp','file'),
  display('Using linprog solver');
  solver = 'linprog';
end

data_mean = column(data_mean);
data_std  = column(data_std);
prior_mean= column(prior_mean);
prior_std = column(prior_std);

% test whether sizes are consistent (otherwise, an error is thrown)
[data_mean, data_std, R];
[prior_mean, prior_std, R'];


if 0,
  % TEST
  prior_mean = ones(5,1);
  prior_std  = ones(5,1);
  data_mean  = ones(4,1);
  data_std   = ones(4,1);
  R = ones(4,5);
  [x_max_posterior, neg_log_posterior, x_min, x_max, x_cov_posterior, y_cov_posterior] = maximise_bayesian_posterior_L1(data_mean, data_std, prior_mean, prior_std, R)
end

% --------------------------------------------------------
% Initialise

n_pars = length(prior_mean);
n_data = length(data_mean);

P = [diag(1./prior_std); ...
     diag(1./data_std) * R];

q = [prior_mean ./ prior_std; ...
     data_mean  ./ data_std];

n_q = length(q);

% --------------------------------------------------------
% from here on, minimise || P * x -q ||_L1
%
% use augmented variables:
% set   d = P * x - q
% split d = dpos - dneg (where dpos >= 0, dneg >= 0)
% set   z = [ dpos; dneg; x ];
% 
% then: minimise sum(dpos) + sum(dneg)  
%       s.t. dpos >= 0, 
%            dneg >= 0, 
%            dpos-dneg = P x - q

f     = [ ones(2*n_q,1); zeros(n_pars,1) ]; 
Aineq = [-eye(n_q), zeros(n_q,n_q+n_pars); ...
          zeros(n_q,n_q),-eye(n_q),zeros(n_q,n_pars)];
bineq = zeros(2*n_q,1);
Aeq   = [eye(n_q), -eye(n_q), -P];
beq   = -q;
lb    = [zeros(2*n_q,1); -x_threshold*ones(n_pars,1)];
ub    = x_threshold * ones(2*n_q+n_pars,1);

switch solver,
  case 'cplex',
    opt = cplexoptimset('Display','off');
    [z_opt,neg_log_posterior, exitflag]  = cplexlp(f, Aineq, bineq, Aeq, beq,lb,ub,[],opt);
  case 'linprog',
    opt = optimset('Display','off');
    [z_opt,neg_log_posterior, exitflag]  = linprog(f, Aineq, bineq, Aeq, beq,lb,ub,[],opt);
end

if exitflag<=0,
  % try once more with relaxed lower and upper bounds
  lb    = [zeros(2*n_q,1); -10^10*x_threshold*ones(n_pars,1)];
  ub    = 10^10*x_threshold * ones(2*n_q+n_pars,1);
%  switch solver,
%    case 'cplex',
      opt = cplexoptimset('Display','off');
      [z_opt,neg_log_posterior, exitflag]  = cplexlp(f, Aineq, bineq, Aeq, beq,lb,ub,[],opt);
%    case 'linprog',
%      opt = optimset('Display','off');
%      [z_opt,neg_log_posterior, exitflag]  = linprog(f, Aineq, bineq, Aeq, beq,lb,ub,[],opt);
%  end
end
  
if exitflag<=0,
  error(sprintf('Solver %s: No solution found - exitflag %d',solver, exitflag)); 
else
  x_max_posterior    = z_opt(2*n_q+(1:n_pars));
  % CHECK:
  % dpos = z_opt(1:n_q);
  % dneg = z_opt(n_q+(1:n_q));
  % [dpos-dneg,P*x-q]
end

% THIS IS THE FORMULA for neg_log_posterior
% neg_log_posterior = sum(abs([x_max_posterior-prior_mean] ./ prior_std)) + sum(abs([R * x_max_posterior - data_mean] ./ data_std))


% --------------------------------------------------------
% now: minimise and maximise every coordinate in x, while keeping the 
% target value below an allowed threshold = optimal target value + 1/2
%
% Why 1/2? In a gaussian distribution, at one std dev distance from the mean,
% the probability goes down by a factor of exp(-1/2); here we want to find the
% lower and upper x values at which, likewise, the probability goes down by a 
% factor of exp(-1/2)

for it = 1:n_pars,
 
  allowed_value_difference = 1/2;
  
  f     = [zeros(2*n_q+n_pars,1)];
  Aineq = [-eye(n_q),       zeros(n_q,n_q+n_pars); ...
            zeros(n_q,n_q), -eye(n_q), zeros(n_q,n_pars); ...
            ones(1,2*n_q),  zeros(1,n_pars)];
  bineq = [zeros(2*n_q,1);...
           neg_log_posterior + allowed_value_difference];

  % find minimal possible x value
  f(2*n_q+it) = 1;
  switch solver,
    case 'cplex',   [z_opt_min,value,exitflag] = cplexlp(f, Aineq, bineq, Aeq, beq,lb,ub,z_opt,opt);
    case 'linprog', [z_opt_min,value,exitflag] = linprog(f, Aineq, bineq, Aeq, beq,lb,ub,z_opt,opt);
  end
  if exitflag<=0, 
    warning(sprintf('Solver %s: No solution found - exitflag %d',solver, exitflag));
    x_min(it,1) = z_opt(2*n_q + it);
  else
    x_min(it,1) = z_opt_min(2*n_q + it);
  end

  % find maximal possible x value
  f(2*n_q+it) = -1;
  switch solver,
    case 'cplex',   [z_opt_max,value,exitflag] = cplexlp(f, Aineq, bineq, Aeq, beq,lb,ub,z_opt,opt);
    case 'linprog', [z_opt_max,value,exitflag] = linprog(f, Aineq, bineq, Aeq, beq,lb,ub,z_opt,opt);
  end
  if exitflag<=0, 
    warning(sprintf('Solver %s: No solution found - exitflag %d',solver, exitflag)); 
    x_max(it,1) = z_opt(2*n_q + it);
  else
    x_max(it,1) = z_opt_max(2*n_q + it);
  end

end

x_cov_posterior = diag([0.5 * [x_max-x_min]].^2);
y_cov_posterior = R*x_cov_posterior*R';
