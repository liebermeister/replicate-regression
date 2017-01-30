function sigma = replicate_regression_insert_stddev(x, sigma, foptions, r)

% sigma = replicate_regression_insert_stddev(x,sigma,foptions)
%
% compute standard deviations from error model

switch foptions.data_scale,
  
  case {'log2', 'log', 'ln', 'log10', 'log2 ratio'},
    %% sigma as absolute standard error of logarithm

    % to find outliers, compare the individual data values to the gene's median for the respective replicate
    for it = unique(r),
      ind_r              = find(r==it);
      x_medians(:,ind_r) = repmat(nanmedian(x(:,ind_r)')',1,length(ind_r));
    end
    % alternative: to find outlier, compare data values to the gene medians (for all replicates together)
    % x_medians          = repmat(nanmedian(x')',1,size(x,2));
    ind_outlier        = find(abs(x - x_medians) > foptions.log_data_adjust_std_threshold);
    sigma              = foptions.data_std_log * ones(size(x));
    sigma(ind_outlier) = foptions.log_data_adjust_std_factor .* abs(x(ind_outlier) - x_medians(ind_outlier));
    
  case 'absolute',
    %% sigma as relative standard error
    sigma = foptions.data_std_relative * x;
    sigma(sigma<foptions.data_std_minimal) = foptions.data_std_minimal;

    if strcmp(foptions.data_scale,'absolute')
      sigma(x>foptions.abs_data_adjust_std_upper) = 2 * sigma(x>foptions.abs_data_adjust_std_upper);
      sigma(x<foptions.abs_data_adjust_std_lower) = 2 * sigma(x<foptions.abs_data_adjust_std_lower);
    end

    
  otherwise
    error('unknown data scale');

end 
