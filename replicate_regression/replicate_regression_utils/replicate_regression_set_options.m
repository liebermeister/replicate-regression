function options = replicate_regression_set_options(t, y, sigma, fix_parameters, options_default, varargin)

% options = replicate_regression_set_options(t, y, sigma, fix_parameters, options_default, varargin)

if length(varargin{1})==1,
  options = varargin{1}{1};
else,
 if mod(length(varargin),2) ==1, 
  options = varargin{1}{end};
 else,
  options = struct;
 end
 for it = 1:floor(length(varargin)/2),
   options = setfield(options,varargin{2*it-1},varargin{2*it});  
 end
end

if ~fix_parameters, 

%% ----------------------------------------------------
%% set options (in structure 'options')

  
options = join_struct(options_default,options);
  
%% --- insert dummy values for missing standard deviations
  
if find(isnan(sigma(find(isfinite(y))))), 
  warning('Standard deviations missing.');
  std_insert = nanmean(column(sigma));
  if isfinite(std_insert), 
    if options.verbose, 
     display('Replacing them by mean standard deviation');  
    end
  else, 
    std_insert = options.insert_std;
    if options.verbose, 
     display(sprintf('Replacing them by values of %f',std_insert)); 
    end
  end
    sigma(isnan(sigma)) = std_insert;
end

if isfinite(options.set_std)
  sigma = options.insert_std * ones(size(sigma));
  if options.verbose, 
    display(sprintf('Replacing all data standard deviation by values of %f',options.set_std)); 
  end
end

%% ----------------------------------------------------
%% make options complete and consistent
%% --------------------------------------------------------


%% --- logaritmic values
  
if options.is_logarithmic==1, options.convert_to_logarithm = 0; end

if options.convert_to_logarithm ==1,
  if length(find(y<0)) + (options.start_value<=0),
    warning('Negative value encountered. Transformation to logarithms is not possible'); 
    options.convert_to_logarithm = 0;
  end
end


%% --- kind of basis functions 

if isfinite(options.start_value),
  options.shift_data           = 'fixed_start_value';
  options.shift_value          = options.start_value;
  options.use_offset           = 0;
  options.basis                = 'sin_half';
  if options.verbose, 
    display('Using fixed start value and sin(pi/2 t/T) basis functions'); 
  end
end

%% --- check and correct the number of components

if length(options.central_mode_width), 
  if size(options.central_mode_width) ~= size(options.deviation_mode_width), 
    error('Prior width vectors have inconsistent lengths');
  else,
    options.n_comp = length(options.central_mode_width);
  end
  if strcmp(options.basis,'cos+sin'), options.n_comp = ceil(options.n_comp/2); end
end

%% --- consider minimal / maximal number of components
  
if options.n_comp < options.n_comp_min, options.n_comp_min = options.n_comp; end
if options.n_comp > options.n_comp_max, options.n_comp_max = options.n_comp; end
  
if ~length(options.central_mode_width),
  
    %% decreasing according to temporal relaxation (see article)
    switch options.basis, 
      case {'polynomial','exp'}, error('Prior widths cannot be computed automatically for this set of basis functions'); 
    end       
    
    T     = max(t);
    if isnan(options.t_smooth), options.t_smooth = T/2; end
    kappa = 1/options.t_smooth;

    switch options.basis, 
      case 'cos',            omega = 2   * pi * [1:options.n_comp_max]       ./ T;
      case 'sin',            omega = 2   * pi * [1:options.n_comp_max]       ./ T;
      case 'cos+sin',        omega = 2   * pi * [1:options.n_comp_max]       ./ T;
      case 'sin_half',       omega = 1/2 * pi * [1:options.n_comp_max]       ./ T;
      case 'sin_horizontal', omega =       pi * [[1:options.n_comp_max]-1/2] ./ T;
    end   

    %% assume two production steps between white noise and protein curves
    mseries = 1./[kappa^2+omega.^2]; 
    mseries = mseries/mseries(1);
    if isnan(options.n_comp), options.n_comp = length(find(mseries>0.05)); 
      if strcmp(options.basis,'cos+sin'), options.n_comp = ceil(options.n_comp/2); end
    end

    %% mseries(1:options.n_comp)

    options.central_mode_width   = options.central_first_mode_width   * mseries(1:options.n_comp);
    options.deviation_mode_width = options.deviation_first_mode_width * mseries(1:options.n_comp);
    
    %% alternative: linearly decreasing mode prior withs
    %% options.central_mode_width   = options.central_first_mode_width   * [options.n_comp:-1:1] / options.n_comp;
    %% options.deviation_mode_width = options.deviation_first_mode_width * [options.n_comp:-1:1] / options.n_comp;
    
    switch options.basis,
      case 'cos+sin',
        options.central_mode_width   = reshape([options.central_mode_width;options.central_mode_width],    1,2*options.n_comp);
        options.deviation_mode_width = reshape([options.deviation_mode_width;options.deviation_mode_width],1,2*options.n_comp);
    end
    
    if isempty(options.central_mode_mean),
      options.central_mode_mean    =  zeros(size(options.central_mode_width));
      options.deviation_mode_mean  =  zeros(size(options.deviation_mode_width));
    end
    
  end
  
  if isfinite(options.t_jump),
    if ~isfinite(options.central_jump_width),
      options.central_jump_width = max(options.central_mode_width);
      options.central_jump_mean  = 0;
    end
    if ~isfinite(options.deviation_jump_width),
      options.deviation_jump_width = max(options.deviation_mode_width);
      options.deviation_jump_mean  = 0;
    end
  end

  end