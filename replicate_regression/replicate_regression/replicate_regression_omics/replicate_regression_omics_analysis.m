% Script replicate_regression_omics_analysis
%
% Run replicate regression analysis for an omics dataset
% and save all results and graphics to files
%
% Input: 
% base_DIR      (directory path for data set, to be referred to in options file]
% foptions_file (filename of options file containing information 
%                       about all further filenames, parameters, etc.)
%
% This script can also be run in a shell:
%
% matlab -r "cd [DIRECTORY]; foptions_file = '[OPTIONS_FILENAME]'; replicate_regression_omics_analysis; quit;"
%
% [DIRECTORY] is the directory in which this script is located
% The omics data filename is given in the options file

eval(default('base_DIR','[]'));

% -------------------------------------------------------------------
% load options

foptions_default = struct(...
'convenience_name', [], ...
    'translation_table_file', [replicate_regression_BASEDIR '/resources/SubtiWiki_id_mapping_20090701_renamed.txt'], ...
    'data_scale',[], ...
    'convert_to_logarithm',[],...
    'log_transformation','arithmetic',...
    'data_time_unit','min',...
    'data_min_data_points',3,...
    'fixed_prior',0,...
    'prior_updating',5,...
    'update_prior_means',0,...
    'run_crossvalidation',1,...
    'postprocess_normalise',0,...
    'graphics_format','png',...
    'graphics_scale','linear',...
    'data_file_csv',[],...
    'result_file_matlab',[],...
    't_smooth',[],...
    'regression_t_interp',[],...
    'regression_tmin',[],...
    'regression_tmax',[],...
    'options_start_value',[],...
    'options_start_at_t',[],...
    'options_constant_before_start',[],...
    'updating_factor',[],...
    'options_out_csv',[],...
    'result_file_csv',[],...
    'graphics_file',[],...
    'log_file',[],...
    'graphics_individual',[],...
    'result_file_zip',[],...
    'normalise_by_median',1,...
    'ignore_std_deviations',1,...
    'data_min_num_replicates',1,...
    'data_std_log',0.17,...
    'data_outliers_threshold',1.5,...
    'data_std_log_outlier',0.7,...
    'data_std_minimal',10^-5,...
    'mark_outliers_percentage',0.02);

% default values for data_std_log, data_outliers_threshold, data_std_log_outlier are for proteome data on log2 scale

foptions = load_options_table(foptions_file, base_DIR);
foptions = join_struct(foptions_default,foptions);

% -------------------------------------------------------------------
% Translate data tables (.csv files) into matlab structures and save them as a .mat file
%
% If necessary, replace (missing/all) standard deviations 
% sigma: Relative measurement error width (according to tests with BaSysBio data)
%        value 0.3: personal communication J. Muntel
%        my own calculation 0.18 (from duplicate data)

cd(foptions.data_dir);

if length(foptions.data_file_csv),
  [data_replicates, data, data_pointwise_average, explanatory_variable] = replicate_regression_load_data_table(foptions.data_file_csv, foptions);
  %% only data_replicates is used for the regression
  
  save(foptions.data_file_matlab, 'data', 'data_replicates', 'data_pointwise_average','explanatory_variable');
  display(sprintf('Saved data to file %s', foptions.data_file_matlab));
end


% ----------------------------------------------------------------
% Run replicate regression (entire omics data set)
% use m-file 'biotable_combine_replicates' (data provided in 'biotable' data format)

if length(foptions.result_file_matlab),

  cd(foptions.data_dir);  
  load(foptions.data_file_matlab);
  options = load_options_table(foptions.options_file);

  if length(foptions.regression_t_interp),
    options.t_interp = foptions.regression_t_interp;
  elseif length(foptions.regression_tmax),
    if ~length(foptions.regression_tmin), foptions.tmin = 0; end
    options.t_interp = foptions.regression_tmin + [foptions.regression_tmax - foptions.regression_tmin] *[0:0.05:1];
  else,
    options.t_interp = [];
  end  
  
  options.run_crossvalidation = foptions.run_crossvalidation; 

  if length(foptions.options_start_at_t),  options.start_at_t  = foptions.options_start_at_t;  end 
  if length(foptions.options_start_value), options.start_value = foptions.options_start_value; end 
  if length(foptions.updating_factor),     options.updating_factor = foptions.updating_factor; end 
  if length(foptions.prior_updating),      options.prior_updating  = foptions.prior_updating;  end
  if length(foptions.t_smooth),            options.t_smooth        = foptions.t_smooth;        end
  if foptions.fixed_prior,                 options.prior_updating = 0; end
  if length(foptions.options_constant_before_start), 
    options.constant_before_start = foptions.options_constant_before_start; 
  end 

  if foptions.postprocess_normalise * [1-strcmp(options.basis,'sin_half')], 
    warning('Normalisation is not possible'); 
  end
  if length(foptions.convert_to_logarithm),
    options.convert_to_logarithm = foptions.convert_to_logarithm; 
  end
  if [foptions.postprocess_normalise * strcmp(foptions.data_scale,'absolute')], 
    options.convert_to_logarithm = 1; 
    display('Switching to log scale');
  end

  options.set_std =   foptions.set_std;
  
  [data_reg, options_complete, options_update] = replicate_regression_complete(data_replicates,options,foptions);

  if foptions.run_crossvalidation,

  crossvalidation_error_mcr_replicate = nanmean([data_reg.combined.DataMean(:) - data_reg.crossvalidation_replicate.DataMean(:)].^2);
  crossvalidation_error_mcr_average   = nanmean([data_reg.combined.DataMean(:) - data_reg.crossvalidation_average.DataMean(:)].^2);
  crossvalidation_error_pooled        = nanmean([data_reg.combined.DataMean(:) - data_reg.pooled_crossvalidation.DataMean(:)].^2);
  crossvalidation_error_naive         = nanmean([data_reg.combined.DataMean(:) - data_reg.naive_crossvalidation.DataMean(:)].^2);
  else
    crossvalidation_error_mcr_replicate = [];
    crossvalidation_error_mcr_average=[];
    crossvalidation_error_pooled = [];
    crossvalidation_error_naive = [];
  end
  
  cd(foptions.result_dir);

  save(foptions.result_file_matlab, 'data_reg', 'options', 'options_complete', 'options_update', 'crossvalidation_error_mcr_average','crossvalidation_error_mcr_replicate', 'crossvalidation_error_pooled', 'crossvalidation_error_naive');
  display(sprintf('Saved regression results to file %s', foptions.result_file_matlab));

end


% ------------------------------------------------
% Show crossvalidation errors

if length(crossvalidation_error_mcr_replicate),
  crossvalidation_error_mcr_replicate
  crossvalidation_error_mcr_average
  crossvalidation_error_pooled
  crossvalidation_error_naive
end

% ------------------------------------------------
% Write results to table files

if length(foptions.options_out_csv),
  cd(foptions.result_dir); 
  save_options_table(options_complete, [foptions.options_out_csv ]);
  display(sprintf('Saving updated options to file %s', [ foptions.options_out_csv ]));
end

display(sprintf('Result file directory %s',foptions.result_dir));

if length(foptions.result_file_csv),
  if strcmp(foptions.result_file_csv(end-3:end),'.tsv'), foptions.result_file_csv = foptions.result_file_csv(1:end-4); end
  if strcmp(foptions.result_file_csv(end-3:end),'.csv'), foptions.result_file_csv = foptions.result_file_csv(1:end-4); end
  cd(foptions.result_dir);
  load(foptions.result_file_matlab);
  replicate_regression_save_data_table(data_reg.average,'R1_R2_R3',[ foptions.result_file_csv '_average.csv' ]);
  fn = fieldnames(data_reg.replicates);
  for it=1:length(fn),
    replicate_regression_save_data_table(data_reg.replicates.(fn{it}),fn{it},[ foptions.result_file_csv '_replicate_' fn{it} '.csv']);
  end
end


% ------------------------------------------------
% Store all results in zip file

if length(foptions.result_file_zip), 
  display('Generating zip file');
  eval(sprintf('! cd %s; cd ..; rm %s',          foptions.data_dir, foptions.result_file_zip));
  eval(sprintf('! cd %s; cd ..; zip %s %s',      foptions.data_dir, foptions.result_file_zip, foptions.data_dir));
  eval(sprintf('! cd %s; cd ..; zip %s %s/*csv', foptions.data_dir, foptions.result_file_zip, foptions.data_dir));
  eval(sprintf('! cd %s; cd ..; zip %s %s',      foptions.data_dir, foptions.result_file_zip, foptions.result_dir));
  eval(sprintf('! cd %s; cd ..; zip %s %s/*csv', foptions.data_dir, foptions.result_file_zip, foptions.result_dir));
  eval(sprintf('! cd %s; cd ..; zip -r %s %s',   foptions.data_dir, foptions.result_file_zip, foptions.graphics_dir));
end


% ------------------------------------------------
% Make graphics and save them

if length(foptions.graphics_file),
  display(sprintf('Graphics directory %s',foptions.graphics_dir));
  cd(foptions.result_dir);
  load(foptions.result_file_matlab);
  cd(foptions.graphics_dir);
  replicate_regression_display_statistics(data_reg, options_complete, foptions.graphics_file, foptions.graphics_format);
end

if length(foptions.graphics_individual),
  cd(foptions.result_dir);
  load(foptions.result_file_matlab);
  cd(foptions.graphics_dir);
  n_elements = size(data_reg.average.DataMean,1);
  ind_show   = 1:n_elements;
  clear gp; 
  gp.convenience_name   = foptions.convenience_name;
  gp.subplot            = [5,5]; 
  gp.log_transformation = foptions.log_transformation; 
  gp.replicate_names    = fieldnames(data_reg.replicates); 
  gp.fontsize           = 10;
  gp.mark_data          = data_reg.presumable_outliers;
  if strcmp(foptions.data_scale,'absolute'), 
    gp.logarithmic_data = 0; 
  else,
    gp.logarithmic_data = 1; 
  end 
  if strcmp(foptions.graphics_scale,'log2'),  
    gp.show_log2 = 1;
  else,
    gp.show_log2 = 0;
  end; 
  gp.image_format = foptions.graphics_format;
  ca;
  biotable_interpolation_graphics_std(data_reg.combined, data_reg.average, data_reg.replicates, gp, ind_show,foptions.graphics_individual);
  if foptions.postprocess_normalise,
    ca;
    biotable_interpolation_graphics_std(data_reg.adjusted.combined, data_reg.adjusted.average, data_reg.adjusted.replicates, gp, ind_show,[foptions.graphics_individual '_adjusted']);
  end
end

