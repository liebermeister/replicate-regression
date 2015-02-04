function base_directory = replicate_regression_DIR()

% This function sets the path of your 'replicate_regression' directory

% Please update the following path upon installation:

base_directory = '/home/wolfram/matlab/wolf_packages/replicate-regression/replicate_regression/';

try
  cd([ base_directory '/resources/']);
catch
  error('Please set the correct path to your replicate_regression directory. You can do this by editing the m-file "replicate_regression_DIR.m"');
end
