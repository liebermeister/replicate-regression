function [data_reg, options_complete, options_update] = replicate_regression_complete(data_replicates,options,foptions)

% [data_reg, options_complete, options_update] = replicate_regression_complete(data_replicates,options)
%
% Multicurve regression for omics data with multiple replicates, including prior updating
%
% 'data_replicates' contains the replicate data sets as fields (all in biotable format, same form!!)
% The entries of the first field are used as unique identifiers
% The data entries are reordered in the combined data set.
%
% data_reg.average, data_reg.replicates, 
% data_reg.central (struct arrays):   regression results
% statistics:                     additional information about the regression procedure
% data_reg.combined (struct arrays) : data structure from intermediate processing steps
%                                 field 'SampleName' contains replicate labels (as numbers)
% data_reg.fit                    fitted data points (reordered!)
% options:                        options for multi-curve regression (see script replicate_regression)
%   with optional additional fields .prior_updating and .updating_factor

% data_collected (struct arrays): data structure from intermediate processing steps

% -----------------------------------------------------------------

eval(default('options','struct','foptions','struct'));
foptions_default = struct('postprocess_normalise',0,'run_crossvalidation',1,'update_prior_means',0);
foptions         = join_struct(foptions_default,foptions);

fn = fieldnames(data_replicates);

data = {};

for it = 1:length(fn),
  data{it} = data_replicates.(fn{it});
end

[data_combined, data_collected] = biotable_join_replicates(data, options);


% -------------------------------------------------------------------------
% If option 'prior_updating' is set, run analysis prior_updating times and update the prior
% Revisit prior assumptions: update data error bars and priors and rerun estimation

run_updating = 0;
if isfield(options,'prior_updating'),
  if isfinite(options.prior_updating),
    if options.prior_updating>0,
      run_updating = 1;
    end
  end
end

time_for_prior_updating_total = nan; 
time_for_prior_updating_per_loop = nan; 

if run_updating,
  my_options = options;  
  tic;
  for it = 1:options.prior_updating,
    display(sprintf('Updating the priors: iteration %d\n',it));
    my_options.verbose             = 0;
    my_options.prior_updating      = nan;
    my_options.run_crossvalidation = 0;
    [data_average, data_rep, data_central, data_cross_average, data_cross_rep, data_fit, options_complete, statistics, options_update] = replicate_regression_biotable(data_combined, my_options);
    my_options = join_struct(my_options, options_update);
  end    
  time_for_prior_updating_total = toc;
  time_for_prior_updating_per_loop = time_for_prior_updating_total/options.prior_updating;  
else,
  options_update = struct;  
end
  
  
% -------------------------------
% Do final replicate regression: build tables data_central and data_average
% and list of tables data_replicates, 

my_options = join_struct(options, options_update);

tic
  display(sprintf('Running the final regression\n'));
  [data_average, data_rep, data_central, data_cross_average, data_cross_rep, data_fit, options_complete, statistics, options_update] = replicate_regression_biotable(data_combined, my_options);
  
  time_for_final_run = toc;

if options.run_crossvalidation,
  data_pooled            = data_combined;
  data_pooled.SampleName = repmat(data_pooled.SampleName(1),length(data_pooled.SampleName),1);

  display(sprintf('Running regression with pooled data for crossvalidation\n'));
  [data_pooled_average, data_pooled_rep, data_pooled_central, data_pooled_crossvalidation ] = replicate_regression_biotable(data_pooled, options);
  options_naive = options;
  
  options_naive.central_offset_width       = inf;
  options_naive.central_first_mode_width   = inf;
  options_naive.central_jump_width         = inf;
  options_naive.deviation_offset_width     = inf;
  options_naive.deviation_first_mode_width = inf;
  options_naive.deviation_jump_width       = inf;

  display(sprintf('Running regression with pooled data and flat prior for crossvalidation\n'));
  [data_naive_average, data_naive_rep, data_naive_central, data_naive_crossvalidation ] = replicate_regression_biotable(data_pooled, options_naive);
end

% -------------------------------------------------------------------------
% prepare output

data_reg.average  = data_average;  
for it = 1:length(fn),
  data_reg.replicates.(fn{it}) = data_rep{it};
end
data_reg.central    = data_central;
data_reg.statistics = statistics;
data_reg.combined   = data_combined;
data_reg.fit        = data_fit;

if options.run_crossvalidation,
  data_reg.crossvalidation_average   = data_cross_average;
  data_reg.crossvalidation_replicate = data_cross_rep;
  data_reg.pooled_crossvalidation    = data_pooled_crossvalidation;
  data_reg.naive_crossvalidation     = data_naive_crossvalidation;
end


% -------------------------------------------------------------------------
% add information about normalised curves and data

if foptions.postprocess_normalise,
  display(sprintf('Running another regression for adjusted curves\n'));
  data_reg.adjusted = replicate_regression_normalise(data_reg, options_complete,foptions.data_scale);
end

% -------------------------------------------------------------------------

if length(foptions.log_file),
  cd(foptions.result_dir);
  file = fopen(foptions.log_file,'w');
  fprintf(file,'Time for prior updating total : %f s\n',time_for_prior_updating_total);
  fprintf(file,'Time for prior updating per loop: %f s\n',time_for_prior_updating_per_loop);
  fprintf(file,'Number of updating loops: %d\n',options.prior_updating);
  fprintf(file,'Number of biological elements: %d\n',size(data_combined.DataMean,1));
  fprintf(file,'Time for final run: %f s\n',time_for_final_run);
  fprintf(file,'Time for final run per element: %f s\n',time_for_final_run/size(data_combined.DataMean,1));
  fclose(file);
end


% ------------------------------------------------
% try to spot outliers based on crossvalidation:
% an upper percentile of absolute mismatches is counted as outliers

% compare, for each data point, the leave-one out fit of the replicate curves 
% to the actual data value and compute the mismatch
% deviation  = abs(data_reg.combined.DataMean-data_reg.crossvalidation_replicate.DataMean);

% compare, for each data point, the leave-one out fit of the replicate curves 
% to the normal (non-leave-one-out) fit and compute the mismatch

if foptions.run_crossvalidation, 
  deviation  = abs(data_reg.crossvalidation_replicate.DataMean - data_reg.fit.DataMean);
  qq         = quantile(deviation(:),1-foptions.mark_outliers_percentage);
  data_reg.presumable_outliers = sparse(deviation>qq);
else
  data_reg.presumable_outliers = [];
end