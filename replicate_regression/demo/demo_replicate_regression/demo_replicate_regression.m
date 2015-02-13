% ----------------------------------------------------------------
% Demo example for Bayesian regression with multiple time series
%
% matlab function used: replicate_regression.m
%
%  1. Create artificial data
%     o True data points (time points 't_true', values 'x_true')
%     o Noisy replicate data points (t,y,l) with labels l for three replicates with systematic errors
%  2. Run replicate_regression
%  3. Display results
%
% Wolfram Liebermeister (2013)
%
% Contact: wolfram.liebermeister@gmail.com

clear; ca


% ---------------------------------
% Create artificial data (three replicates with systematic errors)

[t,y,sigma,r,t_true,x_true, t1,t2,t3,x1,x2,x3,sigma1,sigma2,sigma3,x1_true,x2_true,x3_true,x1_true_all,x2_true_all,x3_true_all] = demo_replicate_regression_create_data;


% ---------------------------------
% Regression options

options                            = struct;
options.convert_to_logarithm       = 0;
options.start_value                = 0.0;
options.std_insert                 = 1.0;
options.central_first_mode_width   = 0.3;
options.deviation_first_mode_width = 0.5;
options.deviation_same_start       = 1;
options.n_comp                     = nan;
options.t_smooth                   = 30;
options.t_jump                     = nan;
options.t_interp                   = -5:1:30;
options.run_crossvalidation        = 1;


% ---------------------------------
% Plot results

fontsize = 24;
colors   = rr_colors;


% ---------------------------------
% True data

figure(1); clf; set(gca,'FontSize',fontsize);

replicate_regression_display(t, y, sigma, r, [],[],[],struct('fignum',1)); hold on;

plot(t_true,x_true,'m--', 'LineWidth',2); hold on
plot(t_true,x1_true_all,'-','Color',colors{1}, 'LineWidth',2);  hold on
plot(t_true,x2_true_all,'-','Color',colors{2}, 'LineWidth',2); hold on
plot(t_true,x3_true_all,'-','Color',colors{3},'LineWidth',2);

axis([0 30 0 1.2]); set(gca,'Fontsize',fontsize); 
xlabel('Time [min]'); ylabel('Protein level [a.u.]'); legend off


% ---------------------------------
% Naive alternative (I) 
% assume that all data points came from one single replicate

options_flat_prior = options;
options_flat_prior.central_first_mode_width   = inf;
options_flat_prior.deviation_first_mode_width = inf;

r_naive      = ones(size(r));

result_naive = replicate_regression(t, y, sigma, r_naive, 0, options);

figure(2); clf; set(gca,'FontSize',fontsize);

replicate_regression_display(t, y, sigma, r, [],[],[],struct('fignum',2)); hold on;
plot(t_true,x_true,'m--', 'LineWidth',2); hold on
plot(result_naive.t,result_naive.x_average,'k-', 'LineWidth',2); hold on
if isfield(result_naive,'x_crossvalidation'),
  plot(t,result_naive.x_crossvalidation,'b*'); hold on
end
axis([0 30 0 1.2]); set(gca,'Fontsize',fontsize); 
xlabel('Time [min]'); ylabel('Protein level [a.u.]'); legend off

result_naive_flat = replicate_regression(t, y, sigma, r_naive, 0, options_flat_prior);

figure(12); clf; set(gca,'FontSize',fontsize);
replicate_regression_display(t, y, sigma, r, [],[],[],struct('fignum',12)); hold on;
plot(t_true,x_true,'m--', 'LineWidth',2); hold on
plot(result_naive_flat.t,result_naive_flat.x_average,'k-', 'LineWidth',2); hold on
if isfield(result_naive_flat,'x_crossvalidation'),
  plot(t,result_naive_flat.x_crossvalidation,'b*'); hold on
end
axis([0 30 0 1.2]); set(gca,'Fontsize',fontsize); 
xlabel('Time [min]'); ylabel('Protein level [a.u.]'); legend off


% ---------------------------------
% Naive alternative (II) 
% Treat replicate curves completely separately and average afterwards

ind1 = (r==1); t1 = t(ind1); y1 = y(ind1); sigma1 = sigma(ind1); r1 = r(ind1);
ind2 = (r==2); t2 = t(ind2); y2 = y(ind2); sigma2 = sigma(ind2); r2 = r(ind2);
ind3 = (r==3); t3 = t(ind3); y3 = y(ind3); sigma3 = sigma(ind3); r3 = r(ind3);

options.central_first_mode_width   = sqrt(options.central_first_mode_width^2+options.deviation_first_mode_width^2);

result_naive1 = replicate_regression(t1, y1, sigma1, r1, 0, options);
result_naive2 = replicate_regression(t2, y2, sigma2, r2, 0, options);
result_naive3 = replicate_regression(t3, y3, sigma3, r3, 0, options);
x_separate_average = [result_naive1.x_average+result_naive2.x_average+result_naive3.x_average]/3;


figure(3); clf; set(gca,'FontSize',fontsize);
h    = replicate_regression_display(t, y, sigma, r, [],[],[],struct('fignum',3)); hold on;
h(4) = plot(t_true,x_true,'m--', 'LineWidth',2); hold on
h(5) = plot(result_naive1.t,result_naive1.x_average,'-', 'Color',colors{1}, 'LineWidth',2);
       plot(result_naive2.t,result_naive2.x_average,'-', 'Color',colors{2}, 'LineWidth',2);
       plot(result_naive3.t,result_naive3.x_average,'-', 'Color',colors{3},'LineWidth',2);
h(6) = plot(result_naive1.t,x_separate_average,     'k-', 'LineWidth',2); 

if isfield(result_naive1,'x_crossvalidation'),
  h(7) = plot(t1,result_naive1.x_crossvalidation,'*', 'Color',colors{1}, 'LineWidth',2);
         plot(t2,result_naive2.x_crossvalidation,'*', 'Color',colors{2}, 'LineWidth',2); 
         plot(t3,result_naive3.x_crossvalidation,'*', 'Color',colors{3},'LineWidth',2); hold off
end

axis([0 30 0 1.2]); set(gca,'Fontsize',fontsize); 
xlabel('Time [min]'); ylabel('Protein level [a.u.]'); legend off


% ---------------------------------
% Replicate regression

[result, options_completed] = replicate_regression(t, y, sigma, r, 0, options);

replicate_regression_display(t, y, sigma, r, [], [], result, struct('fignum',4, 'fontsize',fontsize,'show_crossvalidation',0,'linewidth',2,'show_central',0));

figure(4); axis([0 30 0 1.2]); set(gca,'Fontsize',fontsize); 
xlabel('Time [min]');  ylabel('Protein level [a.u.]'); legend off

replicate_regression_display(t, y, sigma, r, t_true, x_true, result, struct('fignum',5, 'fontsize',fontsize,'show_crossvalidation',0,'linewidth',2,'show_bands',0));

figure(5); axis([0 30 0 1.2]); set(gca,'Fontsize',fontsize); 
xlabel('Time [min]');  ylabel('Protein level [a.u.]'); legend off
