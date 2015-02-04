function sigma = replicate_regression_insert_stddev(x, sigma, foptions, r)

% sigma = replicate_regression_insert_stddev(x,sigma,foptions)
%
% compute standard deviations from error model

switch foptions.data_scale,
  
  case {'log2', 'log', 'ln', 'log10', 'log2 ratio'},
    %% sigma as absolute standard error of logarithm

    % to find outlier, compare data values to the gene medians (for the respective replicate!)
    for it = unique(r),
      ind_r              = find(r==it);
      x_medians(:,ind_r) = repmat(nanmedian(x(:,ind_r)')',1,length(ind_r));
    end
    % alternative: to find outlier, compare data values to the gene medians (for all replicates together)
    % x_medians          = repmat(nanmedian(x')',1,size(x,2));
    ind_outlier        = find(abs(x - x_medians) > foptions.data_outliers_threshold);
    sigma              = foptions.data_std_log * ones(size(x));
    sigma(ind_outlier) = foptions.data_std_log_outlier;
    
  case 'absolute',
    %% sigma as relative standard error
    sigma = foptions.data_std_relative * x;
    sigma(sigma<foptions.data_std_minimal) = foptions.data_std_minimal;
  
  otherwise
    error('unknown data scale');

end 
