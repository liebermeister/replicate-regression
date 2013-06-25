% ----------------------------------------------------------------
% Test example for Bayesian regression with single time series
% (matlab function 'replicate_regression_core_single')
%
% Contents of this script:
%   1. Create artificial data
%     - True data points (time points 't_true', values 'x_true')
%     - Noisy replicate data points (t,x)
%   2. Run replicate_regression
%   3. Display results

clear

% ---------------------------------
% Create artificial data

[t, y, sigma, l, t_true, x_true, t1, t2, t3, x1, x2, x3, sigma1, sigma2, sigma3] = demo_replicate_regression_create_data;


% ---------------------------------
% Run regression

t_new              = t_true;
options            = struct;
options.basis      = 'sin_half';
options.use_offset = 0;
options.n_comp     = 5;
options.mode_width = 0.2 * [2 1 .5 .2 .1];
options.t_smooth   = 2;

[result, parameters, options, sample] = replicate_regression_core_single(t, y, sigma, t_new, options);

% ---------------------------------
% Plot results

replicate_regression_single_display(t, y, sigma, t_true, x_true, result);

figure(1); axis([0 3 0 1.5])
figure(2); axis([0 3 0 1.5])

replicate_regression_single_display([], [], [], [],  [], result.derivative,struct('fignum',[3,4]));

replicate_regression_single_display([], [], [], [],  [], result.production,struct('fignum',[5,6]));
