function all = replicate_regression_save_data_table(data_table, rep_name, basename)

% all = replicate_regression_save_data_table(data_table, rep_name, basename)

eval(default('rep_name','''R1''','basename','[]'));

data_table  = biotable_clean(data_table);
sample_name = data_table.SampleName;
sample_time = data_table.SampleTime;
data_mean   = data_table.DataMean;
data_std    = data_table.DataStd;
info        = data_table.Info;
data_table  = rmfield(data_table,{'DataMean','DataStd','Info','SampleName','SampleTime'});

% -------------------------------------------------

fn  = fieldnames(data_table);
all = {};

[n_items, n_samples] = size(data_mean);

values = reshape([data_mean; data_std],n_items,2*n_samples);
values = roundsd(values, 3);

clear names time type

  for it =1:length(fn),
    all = [all [{['!' fn{it}]}; {''}; {''}; {''}; getfield(data_table,fn{it})]];
  end
  
  all{2,1}='!Time';
  all{3,1}='!Replicate';
  all{4,1}='!ValueType';

  for it=1:n_samples,
    time{1,2*it-1} = num2str(sample_time(it));
    time{1,2*it}   = num2str(sample_time(it));
    if ~isstr(sample_name{it}), sample_name{it} = num2str(sample_name{it}); end 
    names{1,2*it-1}  = [sample_name{it} ' Mean']; 
    names{1,2*it  }  = [sample_name{it} ' Std'];
    types{1,2*it-1}  = 'Mean';
    types{1,2*it  }  = 'Std';
    rep{1,2*it  }    = rep_name;
    rep{1,2*it-1  }  = rep_name;
  end
  all = [all, [names; time; rep; types; num2cell(values)]];

% ------------------------------------------------

if length(basename),
  display(sprintf('Writing file %s', basename));
  table(all,0, basename);
end
