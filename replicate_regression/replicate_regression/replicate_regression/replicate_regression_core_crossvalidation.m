function [result, parameters, options, sample] = replicate_regression_core_crossvalidation(t, y, sigma, r, t_new, flag_single, options);

% [result, parameters, options, sample] = replicate_regression_core_crossvalidation(t, y, sigma, r, t_new, flag_single, options);
%
% Run a crossvalidation loop around the matlab functions 
% 'replicate_regression_core' or 'replicate_regression_core_single'

if flag_single,
  [result, parameters, options, sample] = replicate_regression_core_single(t, y, sigma, t_new, options);
else,
  [result, parameters, options, sample] = replicate_regression_core(t, y, sigma, r, t_new, options);
end

my_options                      = options;
my_options.flag_draw_sample     = 0;
my_options.flag_time_derivative = 0;

result.x_cross_average = nan * y;
result.x_cross_replicate = nan * y;

for it = 1:length(t),
  ind = setdiff(1:length(t),it);
  if flag_single, 
    my_result = replicate_regression_core_single(t(ind), y(ind), sigma(ind), t, my_options);
    result.x_cross_average(it)   = my_result.x(it);
    result.x_cross_replicate(it) = my_result.x(it);
  else,
    my_result = replicate_regression_core(t(ind), y(ind), sigma(ind), r(ind), t, my_options);
    if length(my_result.x_replicate) >= r(it), 
      %% exclude cases where a replicate time series contains only one data point
      result.x_cross_average(it)   = my_result.x_average(it);
      result.x_cross_replicate(it) = my_result.x_replicate{r(it)}(it);
    end
  end
end
