function [result, parameters, options, sample] = replicate_regression_core(t,y,sigma,r,tt,options)

%[x_central, x_replicate, sigma_central, sigma_replicate, parameters, x_average, sigma_average, x_sample] = replicate_regression_core(t,y,sigma,r,tt,options)
%
% Bayesian replicate regression for single time series with multiple replicates
%
% DESCRIPTION
%  Interpolation of multiple time series data from nr replicate experiments
%  The regression functions are represented by linear combinations of basis functions (e.g. Fourier components)
%  The expansion coefficients (=model parameters) given prior distributions and estimated by taking the posterior mode
%
% FUNCTION ARGUMENTS
%  Each data point is a quadruple [t(i), y(i), sigma(i), r(i)] of 
%   - time point t(i) 
%   - measured value y(i) 
%   - standard deviation (error bar) sigma(i)
%   - replicate label r(i) (values from 1,..,nr)
%
%  In this script, the function arguments t, y, sigma, r 
%  are given as ROW vectors (or as matrices, where each row is treated separately)
%
% FUNCTION OUTPUTS 
%   Row vectors of regression curves:
%    'x_central'  central regression curves
%    'x_average'  regression curve averaged over replicates
%    'x_replicate' regression curves for the individual replicates 
%
%   and the corresponding uncertainties (row vectors):
%    'sigma_central', 'sigma_average', and 'sigma_replicate'
%
%  'parameters' is a structure array containing the estimated parameter values 
%        (to be used for statistical evaluation of the prior hyperparameters)
%
%  'x_fit' contains the replicate regression curves, evaluated at the point of original data points
%  'x_sample' is a vector of predicted data, obtained from a random sample from the posterior
%
% OPTIONS FOR THE ALGORITHM
%  Detailed options are given the function argument 'p', 
%  a structure array with (optional) fields: 
%
%    options.basis                 : type of basis functions
%    options.use_offset            : (Boolean) flag for stating that there is a constant basis function
%    options.remove_offset         : (Boolean) flag for stating that the constant basis function should be removed
%    options.n_comp                : number of Fourier components (not including the constant offset)
%                                   (sin and cos for the same wavenumber are counted as one component)
%    options.t_jump                : add (at the beginning) a component that yields a constant offset 
%    options.t_smooth              : time constant for estimating production rates
%    options.constant_before_start : set all basis functions to constant values for negative time values
%
%    options.central_mode_mean    : prior mean for Fourier coefficients alpha for central curve 
%                              (vector; same values are used for for sin and cos)
%    options.central_offset_mean  : the same (scalar) for the constant basis function (default 1)
%    options.central_jump_mean    : the same (scalar) for the jump basis function (default 1)
%    options.deviation_mode_mean  : prior mean for Fourier coefficients beta for deviations from central curve
%    options.deviation_offset_mean: the same, for the constant basis function (default 1)
%    options.deviation_jump_mean  : the same, for the jump basis function (default 1)
%
%    options.central_mode_width    : prior width for Fourier coefficients alpha for central curve 
%                              (vector; same values are used for for sin and cos)
%    options.central_offset_width  : the same (scalar) for the constant basis function (default 1)
%    options.central_jump_width    : the same (scalar) for the jump basis function (default 1)
%    options.deviation_mode_width  : prior width for Fourier coefficients beta for deviations from central curve
%    options.deviation_offset_width: the same, for the constant basis function (default 1)
%    options.deviation_jump_width  : the same, for the jump basis function (default 1)
%
%    options.average_std           : how is the uncertainty of the average curve computed?
%                                    'std_dev_mean':  standard deviation of the mean
%                                    'curve_spread':  standard deviation of the individual curves
%    options.flag_draw_sample      : Draw sample curve parameters and curve from the posterior
%    options.flag_time_derivative : Compute time derivative curves
% 
%
%  The basis functions are adjusted to the final time interval [ta,tb](from tt)
%    'cos'            : cosine function, zero slope at t=ta and t=tb
%    'sin'            : sine function, zero value at t=ta and t=tb
%    'sin_half'       : sine function, zero value at t=ta
%    'sin_horizontal' : sine function, zero value at t=ta, zero slope at t=tb
%    'cos+sin'        : cosine and sine functions, no restriction
%    'polynomial'     : polynomial function, zero value at t=ta
%    'exp'            : exponentially relaxing functions (t<0 => f=0; t>0 => f = 1-exp(t/tau);
%
%  The entire curves are shifted by a constant basis function 
%  This can be suppressed by setting options.use_offset = 0
%
% For an example, see the m-file demo_replicate_regression_core
%
% Wolfram Liebermeister (2010)

% NAMING OF VARIABLES
%   Quantities concerning the prediction (symbols with a hat ^ in the paper) are denoted by 'pred'  

%  ------------------------------------------------------------------------------
% initialisation

if ~length(y),  error('Empty data vector'); end

eval(default('options','struct'));
eval(default('tt','min(t)+[max(t)-min(t)]*[0:0.1:1]'));

% default number of basis functions
if ~isfield(options,'n_comp'), options.n_comp = 3; end 
if strcmp(options.basis,'cos+sin'), 
  options.n_comp = 2*options.n_comp;
end 

options_default = struct(...
    'basis','sin_half', ...
    'use_offset',             0, ...
    'central_offset_mean',    0, ...
    'central_mode_mean',      zeros(options.n_comp,1), ...
    'deviation_offset_mean',  0, ...
    'deviation_mode_mean',    zeros(options.n_comp,1), ...
    'central_offset_width',   1, ...
    'central_mode_width',     ones(options.n_comp,1), ...
    'deviation_offset_width', 1, ...
    'deviation_mode_width',   ones(options.n_comp,1), ...
    't_jump',                 nan, ...
    't_smooth',               nan, ...
    'constant_before_start',  1, ...
    'average_std',            'std_dev_mean', ...
    'flag_draw_sample',       1, ...
    'flag_time_derivative',  0);

options = join_struct(options_default,options);

if isfinite(options.t_smooth), options.flag_time_derivative = 1; end
if length(options.central_mode_width)   ~= options.n_comp, error('wrong number of prior widths'); end
if length(options.deviation_mode_width) ~= options.n_comp, error('wrong number of prior widths'); end
if length(options.central_mode_mean)    ~= options.n_comp, error('wrong number of prior means');  end
if length(options.deviation_mode_mean)  ~= options.n_comp, error('wrong number of prior means');  end

% make sure time vectors are actually rows
t  = column(t)';
tt = column(tt)';

% care for missing sigma values
if isempty(sigma), sigma = ones(size(y)); end

% remove missing values from input vectors
y_orig     = y;
ind_finite = find(isfinite(y));
t          = t(ind_finite);
y          = y(ind_finite);
sigma      = sigma(ind_finite);
r          = r(ind_finite);

nr         = max(r);


% ----------------------------------------------------
% Build matrices of basis functions


[V, V_reg, W, W_reg] = replicate_regression_construct_basis(t, tt, nr, options);


% ----------------------------------------------------
% build complete matrix for all data

% has_label{itt}: bit vector indicating which data points belong to replicate itt

for itt = 1:nr,  has_label{itt} = double(r==itt); end

M         = V;
M_reg     = [];
M_average = V_reg;

for itt = 1:nr,
  M     = [M, diag(has_label{itt}) * V];
  M_reg = matrix_add_block(M_reg,V_reg);
  M_average = [M_average, 1/nr * V_reg];
end

M_reg = [repmat(V_reg,nr+1,1), [zeros(length(tt),size(M_reg,2)); M_reg]];


%  ------------------------------------------------------------------------------
% construct vectors of prior widths
% mu_alpha:    vector of prior means for alpha parameters
% mu_beta:     vector of prior means for beta parameters
% sigma_alpha: vector of prior widths for alpha parameters
% sigma_beta:  vector of prior widths for beta parameters

switch options.basis

  case {'cos','sin','sin_half','sin_horizontal','polynomial','exp'},
    mu_alpha    = [options.central_offset_mean;  column(options.central_mode_mean)];
    sigma_alpha = [options.central_offset_width; column(options.central_mode_width)];
    mu_beta     = repmat([options.deviation_offset_mean;  column(options.deviation_mode_mean)],1,nr);
    sigma_beta  = repmat([options.deviation_offset_width; column(options.deviation_mode_width)],1,nr);

  case 'cos+sin', % duplicate vectors
    mu_alpha = [options.central_offset_mean; ...
         reshape([column(options.central_mode_mean)'; column(options.central_mode_mean)'],2*length(options.central_mode_mean),1)];
    mu_beta = repmat([options.deviation_offset_mean; reshape([column(options.deviation_mode_mean); column(options.deviation_mode_mean)],2*length(options.deviation_mode_mean),1)],1,nr);
    sigma_alpha = [options.central_offset_width; ...
         reshape([column(options.central_mode_width)'; column(options.central_mode_width)'],2*length(options.central_mode_width),1)];
    sigma_beta = repmat([options.deviation_offset_width; reshape([column(options.deviation_mode_width); column(options.deviation_mode_width)],2*length(options.deviation_mode_width),1)],1,nr);
end

if isfinite(options.t_jump),
    mu_alpha    = [mu_alpha;    options.central_jump_mean;];
    sigma_alpha = [sigma_alpha; options.central_jump_width;];
    mu_beta     = [mu_beta;     options.deviation_jump_mean  * ones(1,nr)];
    sigma_beta  = [sigma_beta;  options.deviation_jump_width * ones(1,nr)];
end

mu_beta    = reshape(mu_beta,   prod(size(mu_beta)),1);
sigma_beta = reshape(sigma_beta,prod(size(sigma_beta)),1);


% -----------------------------------------------------
% Bayesian parameter estimation: 
% 1. compute mean and covariance of coefficients
% 2. compute mean and covariance of curves per time point

theta_prior_mean        = [mu_alpha; mu_beta];
theta_prior_cov_inv     = diag(1./[sigma_alpha; sigma_beta].^2);
x_mean                  = y';
x_cov_inv               = diag(1./sigma.^2);

% to avoid problems with ill-conditioned posterior covariance matrix:
epsilon = 10^-5; 

theta_posterior_cov_inv = M' * x_cov_inv * M + theta_prior_cov_inv;
theta_posterior_cov     = inv(theta_posterior_cov_inv + epsilon * eye(size(theta_posterior_cov_inv)) );
theta_posterior_mean    = theta_posterior_cov * [ M' * x_cov_inv * x_mean + theta_prior_cov_inv * theta_prior_mean];

%alternative: (avoiding the matrix inversion needed for theta_posterior_cov)
%theta_posterior_mean    = theta_posterior_cov_inv \ [ M' * x_cov_inv * x_mean];

% -----------------------------------------------------
% curve reconstruction

x_all     = [ M_reg * theta_posterior_mean ]';
sigma_all = sqrt(diag([ M_reg * theta_posterior_cov * M_reg']))';

% central curve and its uncertainty 
x_central     = x_all(1:length(tt));
sigma_central = sigma_all(1:length(tt));

% average curve and its uncertainty
x_average     = [M_average * theta_posterior_mean]';

switch options.average_std,   
  case 'std_dev_mean',
    %% consider the standard deviation of the mean:
    sigma_average    = sqrt( diag( [ M_average * theta_posterior_cov * M_average'] ) )';
  case  'curve_spread',
    %% consider the standard deviation of the individual curves
    sigma_average = std(reshape(x_all(length(tt)+1:end),length(tt),nr)');
end

% individual replicate curves

for itt = 1:nr,
  x_replicate{itt} =  x_all(itt*length(tt)+(1:length(tt)));
  if isempty( x_replicate{itt}),  
    x_replicate{itt} = nan * ones(1,length(tt)); 
    end
  sigma_replicate{itt} =  sigma_all(itt*length(tt)+(1:length(tt)));
  if isempty( sigma_replicate{itt}),  
    sigma_replicate{itt} = nan * ones(1,length(tt)); 
  end
end

% fitted data points

x_fit = [ M * theta_posterior_mean ]';

% -------------------------------------------------------------------------
% model parameters

alpha        = theta_posterior_mean(1:length(sigma_alpha));
beta         = theta_posterior_mean(length(sigma_alpha)+1:end);
beta         = reshape(beta,prod(size(beta))/nr,nr);

parameters.alpha_offset = alpha(1);
parameters.beta_offset  = beta(1,:);
parameters.alpha        = alpha(2:end);
parameters.beta         = beta(2:end,:);
parameters.alpha_jump   = nan;
parameters.beta_jump    = nan;

if isfinite(options.t_jump),
 parameters.alpha        = alpha(2:end-1);
 parameters.beta         = beta(2:end-1,:);
 parameters.alpha_jump   = alpha(end);
 parameters.beta_jump    = beta(end,:);
end

% the following entries possibly refer to the logarithmic values
parameters.y            = y;
parameters.sigma        = sigma;
parameters.residuals    = y - x_fit;

% if necessary, re-insert missing values into x_fit
dummi = x_fit;
x_fit = y_orig;
x_fit(ind_finite) = dummi;

% -------------------------------------------

result.x_central       = x_central;
result.x_average       = x_average;
result.x_replicate     = x_replicate;
result.sigma_central   = sigma_central;
result.sigma_average   = sigma_average;
result.sigma_replicate = sigma_replicate;
result.x_fit           = x_fit;
result.t               = tt;


% --------------------------------------------
% sample one parameter set from posterior and compute the corresponding curves

if options.flag_draw_sample,

theta_sample = theta_posterior_mean + sqrtm(theta_posterior_cov) * randn(size(theta_posterior_mean));
x_sample_all = [M_reg * theta_sample ]';
n_t          = length(result.t);
nr           = max(r);

sample.x_central    = x_sample_all(1:n_t); 
sample.x_replicate  = reshape(x_sample_all(n_t+1:end),n_t,nr);
sample.x_average    = mean(sample.x_replicate,2);
sample.theta        = theta_sample;

alpha        = theta_sample(1:length(sigma_alpha));
beta         = theta_sample(length(sigma_alpha)+1:end);
beta         = reshape(beta,prod(size(beta))/nr,nr);

sample.parameters.alpha_offset = alpha(1);
sample.parameters.beta_offset  = beta(1,:);
sample.parameters.alpha        = alpha(2:end);
sample.parameters.beta         = beta(2:end,:);
sample.parameters.alpha_jump   = nan;
sample.parameters.beta_jump    = nan;

if isfinite(options.t_jump),
  sample.parameters.alpha        = alpha(2:end-1);
  sample.parameters.beta         = beta(2:end-1,:);
  sample.parameters.alpha_jump   = alpha(end);
  sample.parameters.beta_jump    = beta(end,:);
end

end


% --------------------------------------------
% if necessary, compute time derivative and production rate curves

if options.flag_time_derivative
  MW         = W;
  MW_reg     = [];
  MW_average = W_reg;
  for itt = 1:nr,
    MW     = [MW, diag(has_label{itt}) * W];
    MW_reg = matrix_add_block(MW_reg,W_reg);
  MW_average = [MW_average, 1/nr * W_reg];
  end
  MW_reg = [repmat(W_reg,nr+1,1), [zeros(length(tt),size(MW_reg,2)); MW_reg]];
end

if isfinite(options.t_smooth),
  MP         = MW         + 1/options.t_smooth * M        ;
  MP_reg     = MW_reg     + 1/options.t_smooth * M_reg    ;
  MP_average = MW_average + 1/options.t_smooth * M_average;
end

% -----------------------------------------------------
% curve reconstruction (derivatives)

if options.flag_time_derivative

  result.derivative.t = result.t;
result.derivative.x_all                   = [ MW_reg * theta_posterior_mean ]';
result.derivative.sigma_all               = sqrt(diag([ MW_reg * theta_posterior_cov * MW_reg']))';

% central curve and its uncertainty 
result.derivative.x_central     = result.derivative.x_all(1:length(tt));
result.derivative.sigma_central = result.derivative.sigma_all(1:length(tt));

% average curve and its uncertainty
result.derivative.x_average     = [MW_average * theta_posterior_mean]';

switch options.average_std,   
  case 'std_dev_mean',
    %% consider the standard deviation of the mean:
    result.derivative.sigma_average    = sqrt( diag( [ MW_average * theta_posterior_cov * MW_average'] ) )';
  case  'curve_spread',
    %% consider the standard deviation of the individual curves
    result.derivative.sigma_average = std(reshape(result.derivative.x_all(length(tt)+1:end),length(tt),nr)');
end

% individual replicate curves

for itt = 1:nr,
  result.derivative.x_replicate{itt} =  result.derivative.x_all(itt*length(tt)+(1:length(tt)));
  if isempty( result.derivative.x_replicate{itt}),  
    result.derivative.x_replicate{itt} = nan * ones(1,length(tt)); 
    end
  result.derivative.sigma_replicate{itt} =  result.derivative.sigma_all(itt*length(tt)+(1:length(tt)));
  if isempty( result.derivative.sigma_replicate{itt}),  
    result.derivative.sigma_replicate{itt} = nan * ones(1,length(tt)); 
  end
end

end

% -----------------------------------------------------
% curve reconstruction (production rates)

if isfinite(options.t_smooth),

  result.production.t = result.t;
  result.production.x_all                   = [ MP_reg * theta_posterior_mean ]';
  result.production.sigma_all               = sqrt(diag([ MP_reg * theta_posterior_cov * MP_reg']))';
  
% central curve and its uncertainty 
result.production.x_central     = result.production.x_all(1:length(tt));
result.production.sigma_central = result.production.sigma_all(1:length(tt));

% average curve and its uncertainty
result.production.x_average     = [MP_average * theta_posterior_mean]';

switch options.average_std,
  case 'std_dev_mean',
    %% consider the standard deviation of the mean:
    result.production.sigma_average    = sqrt( diag( [ MP_average * theta_posterior_cov * MP_average'] ) )';
  case  'curve_spread',
    %% consider the standard deviation of the individual curves
    result.production.sigma_average = std(reshape(result.production.x_all(length(tt)+1:end),length(tt),nr)');
end

% individual replicate curves

for itt = 1:nr,
  result.production.x_replicate{itt} =  result.production.x_all(itt*length(tt)+(1:length(tt)));
  if isempty( result.production.x_replicate{itt}),  
    result.production.x_replicate{itt} = nan * ones(1,length(tt)); 
    end
  result.production.sigma_replicate{itt} =  result.production.sigma_all(itt*length(tt)+(1:length(tt)));
  if isempty( result.production.sigma_replicate{itt}),  
    result.production.sigma_replicate{itt} = nan * ones(1,length(tt)); 
  end
end

end

parameters.basis_functions.basis = options.basis;
parameters.basis_functions.t     = t;
parameters.basis_functions.V     = V;
parameters.basis_functions.W     = W;
parameters.basis_functions.t_reg = tt;
parameters.basis_functions.V_reg = V_reg;
parameters.basis_functions.W_reg = W_reg;

