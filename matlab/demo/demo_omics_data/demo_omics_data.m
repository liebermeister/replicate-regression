% --------------------------------------------------------------
% This is a test example: a small number of protein measurements 
% from the salt stress experiment (Hahne et al. 2009)
% --------------------------------------------------------------


% --------------------------------------------------------------
% Alternative 1: read all options for omics data analysis from a table file
%                (filename given in variable 'foptions_file')

base_DIR      = [ replicate_regression_DIR 'demo/demo_omics_data/' ];

foptions_file = [ base_DIR '/options/options_demo_omics_data.csv'];

replicate_regression_omics_analysis;


% --------------------------------------------------------------
% Alternative 2: define options in a matlab struct (also called 'foptions_file')
%                The options for the actual replicate regression are again given in a file
%                (whose filename is given by foptions_file.options_file)

base_DIR      = [ replicate_regression_DIR 'demo/demo_omics_data/' ];

foptions_file.data_dir               = [ base_DIR '/data/' ];
foptions_file.result_dir             = [ base_DIR '/results/' ];
foptions_file.graphics_dir           = [ base_DIR '/graphics/' ];
foptions_file.data_file_csv          = 'hahne_salt_stress_TEST_cytosol_data.csv';
foptions_file.data_file_matlab       = 'hahne_salt_stress_TEST_cytosol_data.mat';
foptions_file.options_out_csv        = 'hahne_salt_stress_TEST_cytosol_result_regression_options.csv';
foptions_file.result_file_matlab     = 'hahne_salt_stress_TEST_cytosol_result.mat';
foptions_file.result_file_csv        = 'hahne_salt_stress_TEST_cytosol_result.csv';
foptions_file.graphics_file          = 'hahne_salt_stress_TEST_cytosol_graphics';
foptions_file.log_file               = 'hahne_salt_stress_TEST_cytosol_log.txt';
foptions_file.translation_table_file = [ replicate_regression_DIR '/resources/SubtiWiki_id_mapping_20090701_renamed.txt'];

foptions_file.options_file           = [ replicate_regression_DIR '/resources/mcr_options_nonlogarithmic_data.csv'];
foptions_file.data_time_unit         = 'min';
foptions_file.data_scale             = 'absolute';
foptions_file.normalise_by_median    = 1;
foptions_file.convert_to_logarithm   = 1;
foptions_file.log_transformation     = 'arithmetic';
foptions_file.set_std                = 0.5;
foptions_file.data_std_relative      = 0.25;
foptions_file.data_min_data_points   = 3;

foptions_file.fixed_prior            = 0;
foptions_file.prior_updating         = 3;
foptions_file.updating_factor      = 1.2;
foptions_file.t_smooth               = 30;
foptions_file.update_prior_means     = 0;
foptions_file.run_crossvalidation    = 1;
foptions_file.postprocess_normalise  = 1;

replicate_regression_omics_analysis;
