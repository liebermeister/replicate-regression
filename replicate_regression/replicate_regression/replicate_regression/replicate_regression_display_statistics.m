function replicate_regression_display_statistics( result, options, filename,image_format)

% replicate_regression_display_statistics( result, options, filename,image_format)
%
% Statistical output for replicate time series regression

eval(default('image_format','''eps'''));

%load([data_dir 'salt_stress_data']);
%load([data_dir 'salt_stress_regression_result']);

vv         = result.statistics.alpha_offset.values;
mu         = options.central_offset_mean;
sigma      = options.central_offset_width;
mu_post    = result.statistics.alpha_offset.mean;
sigma_post = result.statistics.alpha_offset.std;
figure(21); histogram_fit(vv, mu, sigma, mu_post, sigma_post); 
set(gca, 'Fontsize',18)
title(sprintf('Distribution of \\alpha_0 values. Prior \\mu: %2.2f, \\sigma: %2.2f', mu, sigma));
xlabel('Constant offset \alpha_0 for central curve'); 

vv         = result.statistics.alpha.values(1,:);
mu         = options.central_mode_mean(1);
sigma      = options.central_first_mode_width;
mu_post    = result.statistics.alpha.mean(1);
sigma_post = result.statistics.alpha.std(1);
figure(22); 
histogram_fit(vv, mu, sigma, mu_post, sigma_post);  
title(sprintf('Distribution of \\alpha_1 values. Prior \\mu: %2.2f, \\sigma: %2.2f', mu, sigma)); 
xlabel('First Fourier mode amplitude \alpha_1 for central curve'); 

vv         = result.statistics.beta_offset.values(:);
vv = vv(vv~=0); % ignore elements from time series with only one replicate
mu         = options.deviation_offset_mean;
sigma      = options.deviation_offset_width;
mu_post    = result.statistics.beta_offset.mean;
sigma_post = result.statistics.beta_offset.std;
figure(23); histogram_fit(vv, mu, sigma, mu_post, sigma_post);  
title(sprintf('Distribution of \\beta_0 values. Prior \\mu: %2.2f, \\sigma: %2.2f', mu, sigma)); 
xlabel('Constant offset \beta_0 for curve deviation'); 

vv  = squeeze(result.statistics.beta.values(1,1,:));
vv = vv(vv~=0); % ignore elements from time series with only one replicate
mu    = options.deviation_mode_mean(1);
sigma = options.deviation_mode_width(1);
mu_post    = result.statistics.beta.mean(1);
sigma_post = result.statistics.beta.std(1);
figure(24); histogram_fit(vv, mu, sigma, mu_post, sigma_post);  
title(sprintf('Distribution of \\beta_1 values. Prior \\mu: %2.2f, \\sigma: %2.2f', mu, sigma)); 
xlabel('First Fourier mode amplitude \beta_1 for curve deviation');

vv    = result.statistics.residuals;
sigma = result.statistics.sigma;
figure(25); clf; plot(sigma,vv,'.','Color',[.6 .6 .6]); hold on
plot([10 0 10],[10 0 -10],'k--'); axis equal; axis([-0.1,2,-2,2]);
rms_relative_residual = sqrt(mean([result.statistics.residuals./result.statistics.sigma].^2));
title(sprintf('Relative residual RMS %f',rms_relative_residual)); 
xlabel('Assumed standard error'); ylabel('Residual');

%% plot residuals
figure(26); clf; 
ind = isfinite(result.combined.DataMean);
q1 = quantile(result.fit.DataMean(ind)-result.combined.DataMean(ind),0.05);
q2 = quantile(result.fit.DataMean(ind)-result.combined.DataMean(ind),0.95);
im(result.fit.DataMean-result.combined.DataMean,[q1,q2]); colorbar
title('Residuals (color scale: 5%-95% quantiles)');

switch image_format
  case 'eps',
    print([ filename '_statistics_Alpha0.eps'], '-f21', '-depsc');
    print([ filename '_statistics_Alpha1.eps'], '-f22', '-depsc');
    print([ filename '_statistics_Beta0.eps' ],  '-f23', '-depsc');
    print([ filename '_statistics_Beta1.eps' ],  '-f24', '-depsc');
    print([ filename '_statistics_Sigma.eps' ],  '-f25', '-depsc');
%    print([ filename '_statistics_Residuals.eps' ],  '-f26', '-depsc');
  case 'png',
    print([ filename '_statistics_Alpha0.png'], '-f21', '-dpng');
    print([ filename '_statistics_Alpha1.png'], '-f22', '-dpng');
    print([ filename '_statistics_Beta0.png' ],  '-f23', '-dpng');
    print([ filename '_statistics_Beta1.png' ],  '-f24', '-dpng');
    print([ filename '_statistics_Sigma.png' ],  '-f25', '-dpng');
%    print([ filename '_statistics_Residuals.eps' ],  '-f26', '-dpng');
end
