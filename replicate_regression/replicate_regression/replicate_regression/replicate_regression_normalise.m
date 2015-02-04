function data_adjusted = replicate_regression_normalise(data_reg, options, data_scale);

% data_adjusted = replicate_regression_normalise(data_reg, options, data_scale);

options.remove_offset = 1;
[reg_average, reg_replicate, reg_central] = replicate_regression_biotable(data_reg.combined,options);

data_adjusted.average    = reg_average;
data_adjusted.replicates = reg_replicate;
data_adjusted.central    = reg_central;

data_adjusted.combined   = data_reg.combined;
fn = fieldnames(data_reg.replicates);
for it = 1:length(fn),
  offsets = data_reg.replicates.(fn{it}).DataMean(:,1);
  ind     = find( cell2mat(data_reg.combined.SampleName) == it);
  switch data_scale,
    case 'absolute',
      data_adjusted.combined.DataMean(:,ind) = data_reg.combined.DataMean(:,ind) ./ repmat(offsets,1,length(ind));
    case {'log2','ln','log','log10','log2 ratio'},
      data_adjusted.combined.DataMean(:,ind) = data_reg.combined.DataMean(:,ind) - repmat(offsets,1,length(ind));
  end
end
