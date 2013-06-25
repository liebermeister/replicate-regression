function [result, parameters, options, sample] = replicate_regression_core_single(t,y,sigma,tt,options)

%[x_central, x_replicate, sigma_central, sigma_replicate, parameters, x_average, sigma_average, x_sample] = replicate_regression_core_single(t,y,sigma,tt,options)
%
% Bayesian replicate regression for a time series with a single replicate
%
% DESCRIPTION
%  Interpolation of single time series data
%  The regression functions are represented by linear combinations of basis functions (e.g. Fourier components)
%  The expansion coefficients (=model parameters) given prior distributions and estimated by taking the posterior mode
%
% FUNCTION ARGUMENTS
%  Each data point is a triple [t(i), y(i), sigma(i)] of 
%   - time point t(i) 
%   - measured value y(i) 
%   - standard deviation (error bar) sigma(i)
%
%  In this script, the function arguments t, y, sigma 
%  are given as ROW vectors (or as matrices, where each row is treated separately)
%
% FUNCTION OUTPUTS 
%   Row vectors of regression curves:
%    'x'  central regression curves
%   and the corresponding uncertainties (row vector):
%    'sigma_x'
%
%  'parameters' is a structure array containing the estimated parameter values 
%        (to be used for statistical evaluation of the prior hyperparameters)
%
%  'x_sample' is a vector of predicted data, obtained from a random sample from the posterior
%
% OPTIONS FOR THE ALGORITHM
%  Detailed options are given the function argument 'p', 
%  a structure array with (optional) fields: 
%
%    options.basis                 : type of basis functions
%    options.use_offset            : (Boolean) flag for stating that there is a constant basis function
%    options.n_comp                : number of Fourier components (not including the constant offset)
%                                    (sin and cos for the same wavenumber are counted as one component)
%    options.t_jump                : add (at the beginning) a component that yields a constant offset 
%    options.t_smooth              : time constant for estimating production rates
%    options.constant_before_start : set all basis functions to constant values for negative time values
%
%    options.mode_mean             : prior mean for Fourier coefficients alpha for regression curve 
%                                      (vector; same values are used for for sin and cos)
%    options.offset_mean           : the same (scalar) for the constant basis function (default 1)
%    options.jump_mean             : the same (scalar) for the jump basis function (default 1)
%    options.mode_width            : prior width for Fourier coefficients alpha for central curve 
%                                      (vector; same values are used for for sin and cos)
%    options.offset_width          : the same (scalar) for the constant basis function (default 1)
%    options.jump_width            : the same (scalar) for the jump basis function (default 1)
%    options.flag_draw_sample      : Draw sample curve parameters and curve from the posterior
%    options.flag_time_derivative  : Compute time derivative curves
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
%  This can be suppressed by setting options.use_offset = 1
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
    'basis',                  'sin_half', ...
    'use_offset',             0, ...
    'offset_mean',            0, ...
    'mode_mean',              zeros(options.n_comp,1), ...
    'offset_width',           1, ...
    'mode_width',             ones(options.n_comp,1), ...
    't_jump',                 nan, ...
    't_smooth',               nan, ...
    'constant_before_start',  1, ...
    'average_std',            'std_dev_mean', ...
    'flag_draw_sample',       1, ...
    'flag_time_derivative',   0);

options = join_struct(options_default,options);

if isfinite(options.t_smooth), options.flag_time_derivative = 1; end
if length(options.mode_width) ~= options.n_comp, error('wrong number of prior widths'); end
if length(options.mode_mean)  ~= options.n_comp, error('wrong number of prior means');  end

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


% ----------------------------------------------------
% Build matrices of basis functions


[M, M_reg, W, W_reg] = replicate_regression_construct_basis(t, tt, 1, options);


%  ------------------------------------------------------------------------------
% construct vectors of prior widths
% mu_alpha:    vector of prior means for alpha parameters
% sigma_alpha: vector of prior widths for alpha parameters

switch options.basis

  case {'cos','sin','sin_half','sin_horizontal','polynomial','exp'},
    mu_alpha    = [options.offset_mean;  column(options.mode_mean)];
    sigma_alpha = [options.offset_width; column(options.mode_width)];

  case 'cos+sin', % duplicate vectors
    mu_alpha = [options.offset_mean; ...
         reshape([column(options.mode_mean)'; column(options.mode_mean)'],2*length(options.mode_mean),1)];
    sigma_alpha = [options.offset_width; ...
         reshape([column(options.mode_width)'; column(options.mode_width)'],2*length(options.mode_width),1)];
end

if isfinite(options.t_jump),
    mu_alpha    = [mu_alpha;    options.jump_mean;];
    sigma_alpha = [sigma_alpha; options.jump_width;];
end


% -----------------------------------------------------
% Bayesian parameter estimation: 
% 1. compute mean and covariance of coefficients
% 2. compute mean and covariance of curves per time point

theta_prior_mean        = [mu_alpha];
theta_prior_cov_inv     = diag(1./[sigma_alpha].^2);
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

x_all                   = [ M_reg * theta_posterior_mean ]';
sigma_all               = sqrt(diag([ M_reg * theta_posterior_cov * M_reg']))';


% central curve and its uncertainty 
x     = x_all(1:length(tt));
sigma = sigma_all(1:length(tt));

% fitted data points
x_fit       = [ M * theta_posterior_mean ]';


% -------------------------------------------------------------------------
% model parameters

alpha        = theta_posterior_mean(1:length(sigma_alpha));

parameters.alpha_offset = alpha(1);
parameters.alpha        = alpha(2:end);
parameters.alpha_jump   = nan;

if isfinite(options.t_jump),
 parameters.alpha        = alpha(2:end-1);
 parameters.alpha_jump   = alpha(end);
end

% the following entries possibly refer to the logarithmic values
parameters.y            = y;
parameters.sigma        = sigma;
parameters.residuals    = y-x_fit;

% if necessary, re-insert missing values into x_fit
dummi = x_fit;
x_fit = y_orig;
x_fit(ind_finite) = dummi;


% -------------------------------------------

result.x       = x;
result.sigma   = sigma;
result.x_fit           = x_fit;
result.t               = tt;


% --------------------------------------------
% sample one parameter set from posterior and compute the corresponding curves

if options.flag_draw_sample,

theta_sample = theta_posterior_mean + sqrtm(theta_posterior_cov) * randn(size(theta_posterior_mean));
x_sample_all = [M_reg * theta_sample ]';
n_t          = length(result.t);

sample.x    = x_sample_all(1:n_t); 
sample.theta        = theta_sample;

alpha        = theta_sample(1:length(sigma_alpha));

sample.parameters.alpha_offset = alpha(1);
sample.parameters.alpha        = alpha(2:end);
sample.parameters.alpha_jump   = nan;

if isfinite(options.t_jump),
  sample.parameters.alpha        = alpha(2:end-1);
  sample.parameters.alpha_jump   = alpha(end);
end

end


% --------------------------------------------
% if necessary, compute time derivative and production rate curves

if options.flag_time_derivative
  MW         = W;
  MW_reg     = W_reg;
end

if isfinite(options.t_smooth),
  MP         = MW         + 1/options.t_smooth * M        ;
  MP_reg     = MW_reg     + 1/options.t_smooth * M_reg    ;
end

% -----------------------------------------------------
% curve reconstruction (derivatives)

if options.flag_time_derivative

  result.derivative.t = result.t;
  result.derivative.x_all                   = [ MW_reg * theta_posterior_mean ]';
  result.derivative.sigma_all               = sqrt(diag([ MW_reg * theta_posterior_cov * MW_reg']))';
  
  %% central curve and its uncertainty 
  result.derivative.x     = result.derivative.x_all(1:length(tt));
  result.derivative.sigma = result.derivative.sigma_all(1:length(tt));

end

% -----------------------------------------------------
% curve reconstruction (production rates)

if isfinite(options.t_smooth),

  result.production.t = result.t;
  result.production.x_all                   = [ MP_reg * theta_posterior_mean ]';
  result.production.sigma_all               = sqrt(diag([ MP_reg * theta_posterior_cov * MP_reg']))';
  
% central curve and its uncertainty 
result.production.x     = result.production.x_all(1:length(tt));
result.production.sigma = result.production.sigma_all(1:length(tt));

end

parameters.basis_functions.basis = options.basis;
parameters.basis_functions.t     = t;
parameters.basis_functions.V     = M;
parameters.basis_functions.W     = W;
parameters.basis_functions.t_reg = tt;
parameters.basis_functions.V_reg = M_reg;
parameters.basis_functions.W_reg = W_reg;
