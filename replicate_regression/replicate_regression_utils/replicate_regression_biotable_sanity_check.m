% [statistics, options_update] = biotable_replicate_regression_sanity_check(options, alpha_offset, alpha, alpha_jump, beta_offset, beta, beta_jump, all_residuals, all_y, all_sigma)
%
% Statistical sanity check for priors in multi-curve regression
% including a suggestion for changing the priors
%
% to check whether spread of posterior mean parameters is consistent 
% with prior spread, iterate this and you get a sort of empirical bayes!

function [statistics, options_update] = biotable_replicate_regression_sanity_check(options, alpha_offset, alpha, alpha_jump, beta_offset, beta, beta_jump, all_residuals, all_y, all_sigma)

flag_jump  = isfinite(options.t_jump);
ind_finite = find(isfinite(all_y));

if strcmp(options.basis,'cos+sin'),
  alpha = reshape(alpha, size(alpha,1)/2,2*size(alpha,2)); 
  beta  = reshape(beta,  size(beta,1)/2,2*size(beta,2),size(beta,3)); 
end

statistics.values            = all_y(ind_finite);
statistics.sigma             = all_sigma(ind_finite);
statistics.residuals         = all_residuals(ind_finite);

statistics.alpha_offset.mean = nanmean(alpha_offset);
statistics.alpha_offset.std  = nanstd(alpha_offset');
statistics.alpha_offset.rms  = sqrt(nanmean(alpha_offset.^2));
statistics.alpha_offset.values = alpha_offset;

statistics.alpha.mean        = nanmean(alpha,2);
statistics.alpha.std         = nanstd(alpha')';
statistics.alpha.rms         = sqrt(nanmean(alpha.^2,2));
statistics.alpha.values      = alpha;

statistics.alpha_jump.mean   = nanmean(alpha_jump);
statistics.alpha_jump.std    = nanstd(alpha_jump');
statistics.alpha_jump.rms    = sqrt(nanmean(alpha_jump.^2));
statistics.alpha_jump.values = alpha_jump;

%% ignore zero values (arising from time series with only 1 replicate)
ind = find(beta_offset ~=0);
statistics.beta_offset.mean  = nanmean(beta_offset(ind));
statistics.beta_offset.std   = nanstd(beta_offset(ind));
statistics.beta_offset.rms   = sqrt(nanmean(beta_offset(ind).^2));
statistics.beta_offset.values= beta_offset;

%% ignore zero values (arising from time series with only 1 replicate)
beta_mask = beta; beta_mask(beta_mask==0) = nan;
statistics.beta.mean         = nanmean(reshape(beta_mask,options.n_comp,prod(size(beta))/options.n_comp),2);
statistics.beta.std          = nanstd( reshape(beta_mask,options.n_comp,prod(size(beta))/options.n_comp),0,2);
statistics.beta.rms          = sqrt(nanmean(reshape(beta_mask.^2,options.n_comp,prod(size(beta))/options.n_comp),2));
statistics.beta.values       = beta;

%% ignore zero values (arising from time series with only 1 replicate)
ind = find(beta_jump ~=0);
statistics.beta_jump.mean    = nanmean(beta_jump(ind));
statistics.beta_jump.std     = nanstd(beta_jump(ind));
statistics.beta_jump.rms     = sqrt(nanmean(beta_jump(ind).^2));
statistics.beta_jump.values  = beta_jump;

display(sprintf('Residuals           (RMS): %f',sqrt(mean(statistics.residuals(:).^2))));
display(sprintf('Standard deviations (RMS): %f',sqrt(mean(statistics.sigma(:).^2))));
display(sprintf('Residual/Std. Dev.  (RMS): %f',sqrt(mean((statistics.residuals(:)./statistics.sigma(:)).^2))));

names = [{'Offset'}; repmat({'Mode'},length(statistics.alpha.std),1);{'Jump'}];

central_prior_vector = [options.central_offset_width;...
			options.central_mode_width';...
			options.central_jump_width];

central_posterior_vector =  [statistics.alpha_offset.std; ...
			     statistics.alpha.std;...
			     statistics.alpha_jump.std];

deviation_prior_vector = [options.deviation_offset_width;...
                          options.deviation_mode_width';...
			  options.deviation_jump_width];

deviation_posterior_vector = [statistics.beta_offset.std; ...
			      statistics.beta.std; ...
			      statistics.beta_jump.std];

  display('Distribution widths for central curves');
  display(print_matrix([central_prior_vector, central_posterior_vector, central_prior_vector<central_posterior_vector],names,{'Prior','Post','Prior too narrow'}))
  
  display('Distribution widths for deviation curves');
  display(print_matrix([deviation_prior_vector, deviation_posterior_vector, deviation_prior_vector<deviation_posterior_vector],names,{'Prior','Post','Prior too narrow'}))
  
if options.verbose,
  
  display(sprintf('To choose the posterior as your new prior, you may set'));
  display(sprintf(' options.central_offset_mean        = %f (old: %f)',statistics.alpha_offset.mean,options.central_offset_mean ));
  display(sprintf(' options.central_offset_width       = %f (old: %f)',statistics.alpha_offset.std,options.central_offset_width));
  display(sprintf(' options.central_first_mode_mean    = %f (old: %f)',statistics.alpha.mean(1),options.central_first_mode_mean ));
  display(sprintf(' options.central_first_mode_width   = %f (old: %f)',statistics.alpha.std(1),options.central_first_mode_width));
  display(sprintf(' options.central_jump_mean          = %f (old: %f)',statistics.alpha_jump.mean,options.central_jump_mean ));
  display(sprintf(' options.central_jump_width         = %f (old: %f)',statistics.alpha_jump.std,options.central_jump_width));
  display(sprintf(' options.deviation_offset_mean      = %f (old: %f)',statistics.beta_offset.mean,options.deviation_offset_mean ));
  display(sprintf(' options.deviation_offset_width     = %f (old: %f)',statistics.beta_offset.std,options.deviation_offset_width));
  display(sprintf(' options.deviation_first_mode_mean  = %f (old: %f)',statistics.beta.mean(1),options.deviation_first_mode_mean ));
  display(sprintf(' options.deviation_first_mode_width = %f (old: %f)',statistics.beta.std(1),options.deviation_first_mode_width));
  display(sprintf(' options.deviation_jump_mean        = %f (old: %f)',statistics.beta_jump.mean,options.deviation_jump_mean ));
  display(sprintf(' options.deviation_jump_width       = %f (old: %f)',statistics.beta_jump.std,options.deviation_jump_width));
  
end

% fit parameters for measurement standard deviation model sigma^2 == a + b * y

if options.update_prior_means,
  options_update.central_mode_mean          = statistics.alpha.mean';
  options_update.central_jump_mean          = statistics.alpha_jump.mean;
  options_update.deviation_mode_mean        = statistics.beta.mean';
  options_update.deviation_jump_mean        = statistics.beta_jump.mean;
else
  options_update.central_mode_mean          = 0 * statistics.alpha.mean';
  options_update.central_jump_mean          = 0 * statistics.alpha_jump.mean;
  options_update.deviation_mode_mean        = 0 * statistics.beta.mean';
  options_update.deviation_jump_mean        = 0 * statistics.beta_jump.mean;
end

options_update.deviation_offset_mean      = statistics.beta_offset.mean;
options_update.central_offset_mean        = statistics.alpha_offset.mean;

options_update.central_offset_width       = options.updating_factor * statistics.alpha_offset.std;
options_update.central_first_mode_width   = options.updating_factor * sqrt(statistics.alpha.std(1)^2 + statistics.alpha.mean(1)^2);
options_update.central_jump_width         = options.updating_factor * statistics.alpha_jump.std;
options_update.deviation_offset_width     = options.updating_factor * statistics.beta_offset.std;
options_update.deviation_first_mode_width = options.updating_factor * sqrt(statistics.beta.std(1)^2 + statistics.beta.mean(1)^2);
options_update.deviation_jump_width       = options.updating_factor * statistics.beta_jump.std;

rms_relative_residual = sqrt(mean((statistics.residuals(:)./statistics.sigma(:)).^2));

if rms_relative_residual < 0.8,
  display(sprintf('Residuals are too small: either error bars are too large or updating factor %f should be decreased\n',options.updating_factor));
end
if rms_relative_residual > 1.2,
  display(sprintf('Residuals are too large: either error bars are too small or updating factor %f should be increased\n',options.updating_factor));
end
