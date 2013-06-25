function [result, options, offsets] = replicate_regression(t, y, sigma, r, flag_fix_parameters, varargin)

% [result, options] = replicate_regression(t, y, sigma, r, varargin)
%
% Bayesian replicate regression for multiple time series measured in replicate
%  - wrapper for the function 'replicate_regression_core'
%  - data are transformed to logarithmic scale if necessary
%
% FUNCTION ARGUMENTS
%
%  t, y, sigma, r:
%    input data (times, values, standard errors, replicate labels)
%    given as row vectors (see replicate_regression_core.m)
%
%  flag_fix_parameters (Boolean, optional):
%    if set to 1, the options given in the following argument(s) will be 
%    accepted without changes (otherwise they will be checked and updated
%
%  varagin (optional) 
%     either list of property/value pairs for algorithm options (list see below).
%     or structure 'options' containing the property/value pairs (this is mandatory if flag_fix_parameters is set to 1)
%     
%     options list is ordered by priority; earlier options override later options
%
%  In converting the data to logarithms, y and sigma can be interpreted either as 
%  *median* and *geometric standard deviation* or as *mean* and *standard deviation* of the data values
%  This is defined by the argument options.transformation
%
% OPTIONS appearing in this function ('X' marks options in function 'replicate_regression_core')
%
%  OPTION                               TYPE     DEFAULT        MEANING
%  options.verbose                      Boolean  1              Output information during regression
%  options.is_logarithmic               Boolean  0              Declare that data are logarithmic
%  options.convert_to_logarithm         Boolean  1              Convert data to logarithms for regression
%  options.log_transformation           string   'arithmetic'   'arithmetic', 'geometric'
%  options.run_crossvalidation          Boolean  0              Run crossvalidation
%  options.set_std                      float    nan            Value to replace all data standard deviations
%  options.insert_std                   float    1              Value to replace missing data standard deviations
%                                                         
%  options.start_at_t                   float    0              Start regression curves at starting time 'start_at_t' (instead of t=0)
%  options.start_value                  float    nan            Fixed start value for regression curves
%  options.shift_data                   string   'mean'         Policy for shifting data before regression
%                                                               {'none', 'fixed_start_value', 'mean', 'initial', 'fixed_1'}
%  options.shift_value                  float    nan            Shift used when shifting the data 
%  options.basis                      X string   'cos+sin'      Type of basis functions
%  options.n_comp                     X int      nan            Fixed number of basis functions
%  options.n_comp_min                   int      1              Minimal number of basis functions
%  options.n_comp_max                   int      20             Maximal number of basis functions
%  options.use_offset                 X Boolean  1              Use constant function as one of the basis functions 
%  options.constant_before_start      X Boolean  0              Set all basis functions constant for t<0
%  options.deviation_same_start         Boolean  0              Enforce identical start values for all replicates
%  options.remove_offset              X Boolean  0              Omit offset when creating the regression curves 
%                                                        
%  options.t_smooth                     float    nan            Time constant for setting decreasing prior widths 
%  options.t_jump                     X float    nan            Time constant for initial jump basis function
%  options.t_interp                     float    t              Time points for interpolated regression curves
%  options.average_std                X string   'std_dev_mean' Type of uncertainty to be reported for average curve
%                                                         
%  options.central_offset_mean        X float    0     Prior mean  sigma_alpha_0    (for alpha_0   )
%  options.central_offset_width       X float    1     Prior width  sigma_alpha_0    (for alpha_0   )
%  options.central_first_mode_mean    X float    0     Prior mean  sigma_alpha_1    (for alpha_1   )
%  options.central_first_mode_width   X float    1     Prior width  sigma_alpha_1    (for alpha_1   )
%  options.central_mode_mean          X vector   []    Prior means sigma_alpha_m    (for alpha_m   )
%  options.central_mode_width         X vector   []    Prior widths sigma_alpha_m    (for alpha_m   )
%  options.central_jump_mean          X float    nan   Prior means sigma_alpha_jump (for alpha_jump)
%  options.central_jump_width         X float    nan   Prior widths sigma_alpha_jump (for alpha_jump)
%  options.deviation_offset_mean      X float    0     Prior mean  sigma_beta_0     (for beta_0    )
%  options.deviation_offset_width     X float    1     Prior width  sigma_beta_0     (for beta_0    )
%  options.deviation_first_mode_mean  X float    0     Prior mean  sigma_beta_1     (for beta_1    )
%  options.deviation_first_mode_width X float    1     Prior width  sigma_beta_1     (for beta_1    )
%  options.deviation_mode_mean        X float    []    Prior means sigma_beta_m     (for beta_m    )
%  options.deviation_mode_width       X float    []    Prior widths sigma_beta_m     (for beta_m    )
%  options.deviation_jump_mean        X float    0     Prior means sigma_beta_jump  (for beta_jump )
%  options.deviation_jump_width       X float    1     Prior widths sigma_beta_jump  (for beta_jump )
%  options.flag_draw_sample           X Boolean  1     Draw sample curve parameters and curve from the posterior
%  options.flag_time_derivative       X Boolean  0     Compute time derivative curves
%
% Wolfram Liebermeister (2011)
%
% wolfram.liebermeister@charite.de

eval(default('flag_fix_parameters','0'));

options_default = struct(...
      'verbose',                    1, ...
      'is_logarithmic',             0, ...
      'run_crossvalidation',        0, ...
      'convert_to_logarithm',       1, ...
      'log_transformation',         'arithmetic', ...
      'set_std',                    nan, ...
      'insert_std',                 1, ...
      'start_at_t',                 0,...
      'start_value',                nan,...
      'shift_data',                 'mean', ...
      'shift_value',                nan,...
      'use_offset',                 1, ...
      'remove_offset',              0, ...
      'central_offset_mean',        0,...
      'central_offset_width',       1,...
      'deviation_same_start',       0, ...
      'deviation_offset_mean',      0,...
      'deviation_offset_width',     1,...
      'basis'        ,              'cos+sin', ...
      'constant_before_start',      0, ...
      'central_mode_mean',          [], ...
      'central_mode_width',         [], ...
      'deviation_mode_mean',        [], ...
      'deviation_mode_width',       [], ...
      'n_comp'   ,                  nan, ...
      'n_comp_min',                 1,...
      'n_comp_max',                 20,...
      't_smooth',                   nan, ...
      't_jump',                     nan, ...
      't_interp'          ,         unique(t), ...
      'average_std',                'std_dev_mean', ...
      'central_first_mode_mean' ,   0 , ...
      'central_first_mode_width' ,  1 , ...
      'deviation_first_mode_mean',  0, ...
      'deviation_first_mode_width', 1, ...
      'central_jump_mean',          0, ...
      'central_jump_width',         1, ...
      'deviation_jump_mean',        0, ...
      'deviation_jump_width',       1 );

options = replicate_regression_set_options(t, y, sigma, flag_fix_parameters, options_default, varargin);

% --------------------------------------------------------
% convert data, run regression, and convert results back
% --------------------------------------------------------


% --- shift time

if options.start_at_t,
  t        = t - options.start_at_t;
  options.t_interp = options.t_interp - options.start_at_t;
end


% --- convert data to logarithms

if options.convert_to_logarithm,
   if options.verbose,  display('Converting data to logarithms'); end
   if isfinite(options.central_offset_width + options.deviation_offset_width + options.central_first_mode_width + options.deviation_first_mode_width + sum(options.central_offset_width)+sum(options.central_mode_width) ),
     if options.verbose, display('I will use the given priors (means and widths) for the logarithmic values'); end
   end
   [y,sigma] = lognormal_normal2log(y,sigma,options.log_transformation);
   options.shift_value = log(options.shift_value);
 end

 
% --- compute offset
% possible options: 'none' 'fixed_start_value' 'mean' 'initial' 'fixed_1'


 switch options.shift_data,
   case 'none',
     options.shift_value = 0;
     options.use_offset  = 1;

   case 'fixed_start_value',
     %% everything has been set above
     options.use_offset  = 0;

   case 'mean',       %% shift by mean value
     options.shift_value = nanmean(y'); 
     options.use_offset  = 1;

   case 'initial',    %% shift by first finite value
     options.shift_value = y(find(isfinite(y))); 
     if options.shift_value, 
       options.shift_value= options.shift_value(1); 
     else 
       options.shift_value=0;
       end
     options.use_offset = 0;

   case 'fixed_1',  
     options.shift_value = 0;  
     options.use_offset = 0;

   otherwise, error(sprintf('Unknown option %s', options.shift_data));
 end

 
% ---

if options.deviation_same_start,
  if options.verbose, display('Same starting point for all replicate curves. Using sin(pi/2 t/T) basis functions'); end
  options.basis                  = 'sin_half';
  options.deviation_offset_width = 10^-5;
end
 
% --- subtract offset
 
y = y - options.shift_value;

% --- run regression

if length(unique(r(find(isfinite(y)))))==1,
  
  %% Workaround for the case that only values for one replicate are
  %% available
  options_single = options;
  options_single.offset_mean = options.central_offset_mean  + options.deviation_offset_mean  ;
  options_single.mode_mean   = options.central_mode_mean    + options.deviation_mode_mean    ;
  options_single.jump_mean   = options.central_jump_mean    + options.deviation_jump_mean    ;
  options_single.offset_width= sqrt(options.central_offset_width.^2 + options.deviation_offset_width.^2);
  options_single.mode_width  = sqrt(options.central_mode_width.^2   + options.deviation_mode_width.^2  );
  options_single.jump_width  = sqrt(options.central_jump_width.^2   + options.deviation_jump_width.^2  );
  if options.run_crossvalidation,
    [result_single, parameters_single, options_single, sample_single] = replicate_regression_core_crossvalidation(t, y, sigma, r, options_single.t_interp, 1, options_single);
  else,
    [result_single, parameters_single, options_single, sample_single] = replicate_regression_core_single(t, y, sigma, options_single.t_interp, options_single);    
    result_single.x_cross_average   = nan * result_single.x;
    result_single.x_cross_replicate = nan * result_single.x;
  end
  sigma(find(~isfinite(y))) = 10 * mean(sigma(find(isfinite(y))));
  y(find(~isfinite(y)))     = nanmean(y);
  if options.run_crossvalidation,
    [result, parameters, options, sample] = replicate_regression_core_crossvalidation(t, y, sigma, r, options.t_interp, 0, options);
  else,
    [result, parameters, options, sample] = replicate_regression_core(t, y, sigma, r, options.t_interp, options);
    result.x_cross_average   = nan * result.x_fit;
    result.x_cross_replicate = nan * result.x_fit;
  end
  result.x_average       = result_single.x;
  result.sigma_average   = result_single.sigma;
  result.x_replicate     = repmat({result_single.x},max(r),1);
  result.sigma_replicate = repmat({result_single.sigma},max(r),1);
  parameters.beta_offset = nan * parameters.beta_offset;
  parameters.beta        = nan * parameters.beta;
  parameters.beta_jump   = nan * parameters.beta_jump;
else,
  if options.run_crossvalidation,
    [result, parameters, options, sample] = replicate_regression_core_crossvalidation(t, y, sigma, r, options.t_interp, 0, options);
  else,
    [result, parameters, options, sample] = replicate_regression_core(t, y, sigma, r, options.t_interp, options);
    result.x_cross_average = nan * result.x_fit;
    result.x_cross_replicate = nan * result.x_fit;
  end
end


% --- add offset

nr = max(r);
x_central          = result.x_central + options.shift_value;
x_average          = result.x_average + options.shift_value;
x_sample_central   = sample.x_central + options.shift_value;
x_sample_average   = sample.x_average + options.shift_value;
x_fit              = result.x_fit + options.shift_value;
x_cross_average  = result.x_cross_average + options.shift_value;
x_cross_replicate  = result.x_cross_replicate + options.shift_value;
x_replicate        = result.x_replicate;
x_sample_replicate = sample.x_replicate;
sigma_central      = result.sigma_central;
sigma_average      = result.sigma_average;
sigma_replicate    = result.sigma_replicate;

for it_r = 1:nr,
  x_replicate{it_r} = x_replicate{it_r} + options.shift_value;
  x_sample_replicate(:,it_r) = x_sample_replicate(:,it_r) + options.shift_value;
end

% --- convert back from logarithm

if options.convert_to_logarithm, 
  [x_central,sigma_central] = lognormal_log2normal(x_central,sigma_central,options.log_transformation);
  [x_average,sigma_average] = lognormal_log2normal(x_average,sigma_average,options.log_transformation);

  x_sample_central  = exp(x_sample_central);
  x_sample_average  = exp(x_sample_average);
  x_fit             = exp(x_fit);
  x_cross_average   = exp(x_cross_average);
  x_cross_replicate = exp(x_cross_replicate);

  for it_r = 1:nr,
   [x_replicate{it_r}, sigma_replicate{it_r}] = lognormal_log2normal(x_replicate{it_r}, sigma_replicate{it_r},options.log_transformation);
   x_sample_replicate(:,it_r) = exp(x_sample_replicate(:,it_r));
  end
  options.shift_value = exp(options.shift_value);
else
  if isfield(result,'derivative'),  derivative = result.derivative; end
  if isfield(result,'production'),  production = result.production; end
end


% --- shift back time

if options.start_at_t,
  t        = t + options.start_at_t;
  options.t_interp = options.t_interp + options.start_at_t;
end


% --- wrap everything into output structure

result                     = struct;
result.t                   = options.t_interp;
result.x_central           = x_central;
result.x_average           = x_average;
result.x_replicate         = x_replicate;
result.sigma_central       = sigma_central;
result.sigma_average       = sigma_average;
result.sigma_replicate     = sigma_replicate;
result.x_sample_central    = x_sample_central;
result.x_sample_average    = x_sample_average;
result.x_sample_replicate  = x_sample_replicate;
result.x_fit               = x_fit;
if options.run_crossvalidation,
  result.x_cross_average = x_cross_average;
  result.x_cross_replicate = x_cross_replicate;
end
result.parameters          = parameters;
result.sample              = sample;

if exist('derivative','var'),  result.derivative = derivative; end
if exist('production','var'),  result.production = production; end

