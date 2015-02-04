function [X_average, X_replicate, X_central, X_cross_average, X_cross_replicate, X_fit, options, statistics, options_update] = replicate_regression_biotable(data, options)

% [X_average, X_replicate, X_central, X_cross_average, X_cross_replicate, X_fit, options, statistics, options_update]
%     = replicate_regression_biotable(data, options)
%
% Multicurve regression for omics data set with a fixed prior
%
% INPUTS
% data:    biotable format
%          data.SampleName must contain number labels of data replicates
% options: structure with directives for replicate regression 
%
% OUTPUTS 
% X_central, X_replicate, X_average (biotables)
%
% Extra entries in options (that are not present in replicate_regression)
% options.set_rel_std
% options.set_abs_std
% options.logarithmic_data
% options.jump_names (list of strings):   names of data items reuireing a 'jump' basis function in multi-curve regression


% ----------------------------------------------------
% initialise

options_default = struct('verbose',0,'update_prior_means',0,'run_crossvalidation',1,'is_logarithmic',nan,'convert_to_logarithm',1,...
		   'shift_data','none','flag_remove_negative',0,...
		   'basis','cos+sin','t_jump',nan,'insert_std',1,'n_comp',5,'set_std',nan,...
                   'start_at_t',min(data.SampleTime));

options_default.jump_names = {};

options = join_struct(options_default, options);

dummi = data;
if isfield(dummi,'logData'), dummi = rmfield(dummi,{'logData','file','metainfo'}); end
dummi = rmfield(dummi,{'DataMean','DataStd','Info','SampleName','SampleTime'});
fn    = fieldnames(dummi);

names = getfield(data,fn{1});

nr = length(unique(cell2mat(data.SampleName)));

% ----------------------------------------------------
% handle time values

if isempty(options.t_interp),
%  options.t_interp = unique(data.SampleTime);
  options.t_interp = min(data.SampleTime) + [max(data.SampleTime)-min(data.SampleTime)] * [0:0.05:1];
end

% ----------------------------------------------------
% handle negative values

if options.flag_remove_negative,
  if find(data.DataMean(data.DataMean<=0)),
  display('Removing negative values');
  data.DataStd(data.DataMean <=0) = nan;
  data.DataMean(data.DataMean<=0) = nan;
  end
end

if options.convert_to_logarithm,
  if length(find(data.DataMean<0)),
    warning('Negative value encountered. Transformation to logarithms is not possible'); 
    options.convert_to_logarithm = 0;
  end
end

if isfinite(options.set_std),
  display(sprintf('Replacing all standard deviations by a fixed value %f',options.set_std));
  data.DataStd = options.set_std * ones(size(data.DataStd));
else,
  
  if find(isnan(data.DataStd(find(isfinite(data.DataMean))))),
    display(sprintf('Standard deviations missing. Replacing them by a predefined value %f', options.insert_std));
    data.DataStd(find(isnan(data.DataStd))) = options.insert_std; 
  end
  
  if find(data.DataStd==0),
    display(sprintf('Vanishing standard deviations found. Replacing them by a predefined value %f', options.insert_std));
    data.DataStd(find(data.DataStd==0)) = options.insert_std; 
  end
  
end


% ----------------------------------------------------
% Build table 'X_central' with regression_central data 

data_empty            = data;
data_empty.SampleName = {}; 
data_empty.SampleTime = []; 
data_empty.DataMean   = [];
data_empty.DataStd    = [];
X_central = data_empty;

for itf = 1:length(fn),  
  X_central = setfield(X_central,fn{itf},getfield(data,fn{itf}) ); 
end
X_central.SampleName = num2cell(options.t_interp');
X_central.SampleTime = options.t_interp';

for it_r = 1:nr,
  X_replicate{it_r} = X_central;
end

X_average         = X_central;

if options.run_crossvalidation, 
  X_cross_average          = data;
  X_cross_average.DataMean = nan * X_cross_average.DataMean;
  X_cross_average.DataStd  = X_cross_average.DataMean;
else,
  X_cross_average = [];
end
X_cross_replicate =  X_cross_average;

alpha_offset  = [];
alpha         = [];
beta_offset   = [];
beta          = [];
all_residuals = [];
all_y         = [];
all_sigma     = [];


% ------------------------------------------------------------------------------
% Call replicate_regression once to fix the parameters (in structure 'options')
% (same for all regression runs)

[result, options] = replicate_regression(data.SampleTime', data.DataMean(1,:), data.DataStd(1,:), cell2mat(data.SampleName)', 0, options);

X_fit = data;
X_fit.DataStd = nan * X_fit.DataStd;
  

for it=1:length(names),
  
  if mod(it,50)==0, 
    display(sprintf('%d/%d',it,length(names))); 
  end

  %% if option 'options.jump_names' is given, then ONLY the indicated elements can show a jump
  
  options.jump_names = options.jump_names(find(cellfun('length',options.jump_names)));
  my_t_jump = options.t_jump;
  if length(options.jump_names), 
    if ~length(find(strcmp(names{it},options.jump_names))), my_t_jump = nan; end
  end
  
  t     = data.SampleTime';
  y     = data.DataMean(it,:);
  sigma = data.DataStd(it,:);
  r     = cell2mat(data.SampleName)';

  %% remove all missing data points and make sure that each replicate contains at least one data point
  
  ind_data_present   = [];
  accept_data_points = [];
  for itt = 1:length(r),
    if sum(isfinite(y(r==itt))), 
      ind_data_present   = [ind_data_present itt];
      accept_data_points = [accept_data_points find(r == itt .* isfinite(y))];
    end
  end
  t     = t(accept_data_points);
  y     = y(accept_data_points);
  sigma = sigma(accept_data_points);
  r     = r(accept_data_points);
  if ~length(t), warning(sprintf('No data point available for item %s', names{it})); end

  result = replicate_regression(t, y, sigma, r, 1, 't_jump', my_t_jump, options);
  
  X_central.DataMean(it,:) = result.x_central;
  X_central.DataStd(it,:)  = result.sigma_central;
  X_average.DataMean(it,:) = result.x_average;
  X_average.DataStd(it,:)  = result.sigma_average;

  if options.run_crossvalidation,
    X_cross_average.DataMean(it,accept_data_points)   = result.x_cross_average;
    X_cross_replicate.DataMean(it,accept_data_points) = result.x_cross_replicate;
  end

  for it_r = 1:nr,
    X_replicate{it_r}.DataMean(it,:) = X_central.DataMean(it,:);
    X_replicate{it_r}.DataStd(it,:)  = X_central.DataStd(it,:);
    if find(it_r == ind_data_present),
      X_replicate{it_r}.DataMean(it,:) = result.x_replicate{find(it_r == ind_data_present)};
      X_replicate{it_r}.DataStd(it,:)  = result.sigma_replicate{find(it_r == ind_data_present)};
      X_fit.DataMean(it,accept_data_points) = result.x_fit;
    else
      X_replicate{it_r}.DataMean(it,:) = nan;
      X_replicate{it_r}.DataStd(it,:)  = nan;
    end    
  end

  %% ---------------------------------------------------
  %% collect values needed for statistical sanity check
  
  %% variant 1: use posterior means of curve parameters
  %% parameters        = result.parameters; 
  
  %% variant 2: use curve parameters sampled from the posterior
  parameters        = result.sample.parameters; 

  all_y             = [all_y,         y];
  all_sigma         = [all_sigma,     sigma];
  all_residuals     = [all_residuals, y - result.x_fit];
  alpha_offset(it)  = parameters.alpha_offset;
  alpha(:,it)       = parameters.alpha;
  beta_offset(it,ind_data_present) = parameters.beta_offset(ind_data_present);
  beta(:,ind_data_present,it) = parameters.beta(:,ind_data_present);
  alpha_jump(it)    = parameters.alpha_jump;
  beta_jump(it,:)   = parameters.beta_jump;

end

% ----------------------------------------------------------------------
% statistical sanity check
%
% to check whether spread of posterior mean parameters is consistent 
% with prior spread, iterate this and you get a sort of empirical bayes!

[statistics, options_update] = replicate_regression_biotable_sanity_check(options, alpha_offset, alpha, alpha_jump, beta_offset, beta, beta_jump, all_residuals, all_y, all_sigma);

statistics.basis_functions = result.parameters.basis_functions;
