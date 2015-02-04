% -----------------------------------------------------------------------------------------------------------
% This is a test example: a small number of protein measurements from the Hahne et al. salt stress experiment
% -----------------------------------------------------------------------------------------------------------

% Note that the string 'replicate_regression_DIR' points to the directory 'replicate_regression' 
% of your installation. You need to set this path by editing the m-file 'replicate_regression_DIR'

% Directory name for the omics set
base_DIR      = [ replicate_regression_DIR '/demo/omics_data_example/' ];

% Name of options file
foptions_file = [ base_DIR '/options/options_omics_data_example.csv'];

% Run script for replicate regression of omics data
replicate_regression_omics_analysis;
