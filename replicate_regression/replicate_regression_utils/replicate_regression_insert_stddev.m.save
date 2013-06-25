function sigma = replicate_regression_insert_stddev(x, sigma, foptions)

% sigma = replicate_regression_insert_stddev(x,sigma,foptions)
%
% compute standard deviations from error model

switch foptions.data_scale,
  case 'log2 ratio',
    %% sigma as absolute standard error of logarithm
    sigma = foptions.data_std_abs * ones(size(x));
    sigma(find(abs(x)>foptions.data_outliers_threshold)) = foptions.data_std_abs_outlier;
    
  case 'absolute',
    %% sigma as relative standard error
    sigma = foptions.data_std_relative * x;
  
  otherwise
    error('unknown data scale');

end 
