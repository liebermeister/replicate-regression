function replicate_regression_omics(data_file, foptions_file, base_directory, goptions_file)

% replicate_regression_omics(data_file, foptions_file, base_directory)
%
% Bayesian replicate regression for omics data
%
% FUNCTION ARGUMENTS:
%
%  data_file:                omics data file (full directory path)
%
%  foptions_file:            table file containing the options (full directory path)
%
%  base_directory:           directory name for results (full directory path)
%
%  goptions_file (optional): table file containing the options for 
%                            single protein graphics (full directory path)
% 
% OUTPUT: Data and graphics are written to files
%
% --------------------------------------------------------------------------------
% Format of data file (tab-separated text, see examples):
%
%   Line 1: Headers
%           Headers of protein names columns (e.g., !BSUnumber, !BGnumber, !GeneName, !UniprotID),
%           followed by sample names (as headers of data columns)
%   Line 2: Time points
%           first column: !Time       data columns: time points (numbers)
%   Line 3: Replicate numbers
%           first column: !Replicate  data columns: replicate names
%   Line 4 (optional): !ValueType ('Value', 'Mean', or 'Std')
%   Further lines: numerical data
%
%
% --------------------------------------------------------------------------------
% Format of options file (tab-separated text, see examples):
%   Each line contains one attribute:
%     first column: attribute name
%     second column: attribute value (string or number)
%     all further columns are ignored
%
%   Lines starting with the '%' character are ignored (can be used for comments)
%   The attribute 'options_file' allows to declare another options file containing default options 
%   The attribute 'data_file_csv' contains the name of the data file
%
%
% --------------------------------------------------------------------------------
% Attributes in options file (for replicate_regression_omics_analysis)
%
%   data_dir                directory name for data files
%   result_dir              directory name for result files
%   graphics_dir            directory name for graphics
%   data_file_csv           filename for data file (tsv format, see examples)
%   data_file_matlab        filename for matlab data file (written during the analysis)
%   options_file            filename for default options file (tsv format)
%   options_out_csv         filename for completed options file (tsv format, written during analysis)
%   translation_table_file  filename for ID mapping table (see example)
%   result_file_matlab      filename for 
%   result_file_csv         filename for        hahne_salt_stress_cytosol_result.tsv
%   result_file_zip         filename for        hahne_salt_stress_result.zip
%   graphics_file           file basename  for graphics
%   data_time_unit          time unit ('min')
%
%   data_scale              'absolute' or 'log2' (also 'ln','log','log10','log2 ratio'; these are all treated like 'log2');
%   data_min_num_replicates minimal number of valid replicates (genes with less valid replicates are discarded; default 1)
%  For "absolute" data:
%    abs_data_adjust_std_upper     upper threshold; points above are outliers (increase std dev by factor of 3)
%    abs_data_adjust_std_lower     lower threshold; points below are outliers (increase std dev by factor of 3)
%    data_std_relative       default for relative standard deviation
%    data_std_minimal        minimal standard deviation
%  For logarithmic data ( 'log2', 'ln','log','log10','log2 ratio')
%    data_std_log                  default for standard deviation (on log scale)
%    log_data_adjust_std_threshold threshold for data values (on chosen log scale) for which std dev is modified
%                                  (criterion:  | [data value] - [median for this gene & replicate] | > threshold )
%                                  for the inserted std dev, see next entry
%    log_data_adjust_std_factor    new std width = factor * absolute deviation from median
%
%   data_min_data_points    minimal number of data points required in the analysis (default 3)
%                           at least one replicate has to reach this number, points are times t<0 do not count
%                           replicates with less data points are ignored
%   convert_to_logarithm    convert (nonlogarithmic) data to logarithms for replicate regression (Boolean)
%   log_transformation      type of transformation 
%                           'arithmetic': data=mean values and plotting on absolute scale
%                           'geometric' : data = median values and plotting on log scale (but data on absolute scale)
%   ignore_std_deviations   Boolean, ignore standard deviations given in data
%   basis                   string: type of basis functions
%   fixed_prior             keeping the prior fixed? (Boolean, default 0)
%   prior_updating          number of prior updating iterations (default 10)
%   updating_factor         prior updating factor, default 1.2
%   updating_factor_final   prior updating factor before last regression; default [] (not set)
%   update_prior_means      change parameter means from 0 to posterior means while updating? default 0
%   t_smooth                time constant defining how prior widths depend on the frequency
%   options_start_value     fixed starting value -> to be inserted into options as options.start_value
%   options_start_at_t      Starting time point for changes (after constant behaviour) to be inserted into options as options.start_at_t
%   options_constant_before_start Boolean (keep curves constant before starting time) to be inserted into options as options.constant_before_start
%   regression_t_interp     time points for regression (optional)
%   regression_tmin         start time for regression (optional)
%   regression_tmax         end time for regression (optional)
%   crossvalidation         run crossvalidation? (Boolean, default 0)
%   postprocess_normalise   Boolean, default 1
%   graphics_individual     file basename (used in script replicate_regression_omics_selected'
%   graphics_scale          default 'log2', 'linear'
%   graphics_format         'eps', 'png' (for technical reasons, 'eps' needs to be written in single quotes)
%   convenience_name        type of protein names to be used in graphics (default 'SubtiWiki_20090701')
%   normalise_by_median     (Boolean, default TRUE)
%   mark_outliers_percentage percentage of data points to be marked as outliers based on crossvalidation error
%
% Additional attributes in options file for individual graphics 
%   (function 'replicate_regression_omics_selected')
%
%   graphics_scale          'log2','linear'
%   postprocess_normalise   1
%   element_id              id (or list, selected by |)
%   element_name            name (or list, selected by |)
%   delimiter_symbol        symbol for delimiting list of elements (in element_id, element_name)
%   title_string            title for graphics
%   x_label                 x label for graphics
%   y_label                 y label for graphics
%   plot_data               produce plot for data (single element)
%   plot_replicates         produce plot with replicates (single element)
%   plot_regression         produce plot for regression curves (single element)
%   plot_all                produce joint plots for all elements

foptions.data_file_csv = data_file;

replicate_regression_omics_analysis;

if exist('goptions_file','var'),
  foptions_file          = goptions_file;
  replicate_regression_omics_selected;
end

ca;