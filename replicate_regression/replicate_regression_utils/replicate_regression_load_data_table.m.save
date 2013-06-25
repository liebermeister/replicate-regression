function [data, replicate_data, data_pointwise_average, explanatory_variable] = replicate_regression_load_data_table(filename, foptions)

%filename = '/home/wolfram/projekte/data_integration/replicate_regression/examples/otto_glucose_starvation/data/otto_BEF_data.tsv';

% ------------------------------------------------------------------
% load data table and create structure 'replicate_data'

foptions_default = struct('data_std_log', 0.3, 'data_std_relative', 0.3, 'translation_table_file',[], 'convenience_name', [], 'data_scale', [], 'data_outliers_upper',inf,'data_outliers_lower',-inf);

foptions_default.ID_potential = {'GiNumber','BSUnumber','UniprotID','GeneName','BGnumber'};

eval(default('foptions','struct'));
foptions = join_struct(foptions_default,foptions);

Data = sbtab_table_load(filename);

fn   = fieldnames(sbtab_table_get_all_columns(Data)); primary_column = fn{1};

% ----------------------------------------------------------------

explanatory_variable = fieldnames(Data.data.attributes);
explanatory_variable = explanatory_variable{1};

if ~strcmp('Time',  explanatory_variable),
  display(sprintf('Explanatory variable %s will be formally treated as time.',explanatory_variable))
end


% ----------------------------------------------------------------
%WARNING: the treatment of means and standard deviations has not been tested yet!

if isfield(Data.data.attributes,'ValueType'),
  ind_sample_value = find(strcmp(Data.data.attributes.ValueType,'Value') + strcmp(Data.data.attributes.ValueType,'Mean') );
  ind_sample_std   = find(strcmp(Data.data.attributes.ValueType,'Std'));

  X                = cell_string2num(Data.data.data);
  Values           = X(:,ind_sample_value);
  StdDev           = nan * Values;
  SampleName       = Data.data.headers(ind_sample_value);
  SampleTime       = cell_string2num(Data.data.attributes.Time(ind_sample_value));
  SampleReplicate  = Data.data.attributes.Replicate(ind_sample_value);
  
  Std_Values          = X(:,ind_sample_std);
  Std_SampleName      = Data.data.headers(ind_sample_std);
  Std_SampleTime      = cell_string2num(Data.data.attributes.Time(ind_sample_std));
  Std_SampleReplicate = Data.data.attributes.Replicate(ind_sample_std);


  for it1 = 1:length(ind_sample_std);
    for it2 =  1:length(ind_sample_value)
      if strcmp(SampleName{it2},Std_SampleName{it1}) * [SampleTime(it2) == Std_SampleTime(it1)] * strcmp(SampleReplicate{it2},Std_SampleReplicate{it1}),
        StdDev(:,it2) = Std_Values(:,it1);
      end
    end
  end
  
else,
  
  Values           = cell_string2num(Data.data.data);
  SampleName       = Data.data.headers;
  SampleTime       = cell_string2num(Data.data.attributes.Time);
  SampleReplicate  = Data.data.attributes.Replicate;
  StdDev           = nan * Values;

end 

% --------------------------------------------------
% filter out all elements in which no replicate reaches the 
% minimal number of data points

all_rep = unique(SampleReplicate);
r       = label_names(SampleReplicate,all_rep)';

dd = zeros(size(Values,1),1);

for it = 1:max(r),
  dd = dd + [sum( isfinite(Values(:,find([r==it] .* [SampleTime>=0] ))) ,2) >= foptions.data_min_data_points];
end

% display('Warning: using only first 100 Proteins');
% data_available = data_available(1:100);

data_available = find(dd);
Values         = Values(data_available,:);
StdDev         = StdDev(data_available,:);

% remove outlier values
%Values(Values>foptions.data_outliers_upper) = nan;
%Values(Values<foptions.data_outliers_lower) = nan;

if strcmp(foptions.data_scale,'absolute') * length(find(Values<0)),
  error('Negative values encountered');
end

StdDev_guess = replicate_regression_insert_stddev(Values,[],foptions);
StdDev(isnan(StdDev)) = StdDev_guess(isnan(StdDev));

% ---------------------------------------------------------------
% normalise by median per sample

switch foptions.data_scale,
  case 'absolute',
    dum = repmat(nanmedian(Values),size(Values,1),1);
    Values = Values ./ dum;
    StdDev = StdDev ./ dum;
  case {'log2','ln','log','log10','log2 ratio'},
    Values = Values - repmat(nanmedian(Values),size(Values,1),1);
  otherwise,
    error('unknown data scale');
end

% ---------------------------------------------------------------

ID_found = fieldnames(sbtab_table_get_all_columns(Data));
if isempty(foptions.convenience_name), foptions.convenience_name = ID_found{1}; end

ll = label_names(foptions.ID_potential,ID_found);
ID_list = foptions.ID_potential(find(ll));

for it = 1:length(ID_list),
  if isfield(sbtab_table_get_all_columns(Data),ID_list{it}), 
    dummi = sbtab_table_get_column(Data,ID_list{it})
    d.(ID_list{it}) = dummi(data_available);  
  end
end

% if the preferred ID is not yet present

if ~isfield(Data.column.column,foptions.convenience_name),
  translation_table = load_any_table(foptions.translation_table_file);
  if strcmp(foptions.convenience_name,'BSUnumber'),
    d.BSUnumber     = UniProt_to_bsu(d.UniprotID,translation_table);
  else,
    d.(foptions.convenience_name) = repmat({''},length(d.(ID_list{1})),1);
    from = label_names(ID_list(1),translation_table(1,:));
    to   = label_names({foptions.convenience_name},translation_table(1,:));
    translation_from   = translation_table(:,from);
    translation_to     = translation_table(:,to);
    ll                = label_names(lower(d.(ID_list{1})),lower(translation_from));
    d.(foptions.convenience_name)(find(ll),1) = translation_to(ll(find(ll)));
  end
end


% ---------------------------------------------------------------

d.SampleName  = {};
d.SampleTime  = [];
d.DataMean    = [];
d.DataStd     = [];
d.Info        = {};%'Data type', 'Protein ratio'; 'Unit', '1'; 'Time unit', 'min'};

datamean = [];
datastd  = [];

for it = 1:length(all_rep),
  ind = find(strcmp(all_rep{it},SampleReplicate));
  replicate_data.(all_rep{it})            = d;
  replicate_data.(all_rep{it}).SampleName = SampleName(ind)';
  replicate_data.(all_rep{it}).SampleTime = SampleTime(ind)';
  replicate_data.(all_rep{it}).DataMean   = Values(:,ind);
  replicate_data.(all_rep{it}).DataStd    = StdDev(:,ind);
  datamean(:,:,it) = Values(:,ind);
  datastd(:,:,it)  = StdDev(:,ind);
end

% attention: this works only if all replicates refer to the same time points
data_pointwise_average = d;
data_pointwise_average.SampleName = replicate_data.(all_rep{1}).SampleName;
data_pointwise_average.SampleTime = replicate_data.(all_rep{1}).SampleTime;
data_pointwise_average.DataMean = squeeze(nanmean(datamean,3));
data_pointwise_average.DataStd  = sqrt(squeeze(nanmean(datastd.^2,3)));


% ---------------------------------------------------------------
% build data structure 'data'
% remove proteins for which no data are available

data = {};

for ind = 1:length(d.(foptions.convenience_name));
  my_id = d.(foptions.convenience_name){ind};
  if length(my_id),
    if length(str2num(my_id)), my_id = ['ID' my_id]; end 
    data.(my_id).(foptions.convenience_name) = my_id;
    for it =1:length(ID_list),
      data.(my_id).(ID_list{it}) = d.(ID_list{it}){ind};
    end
    data.(my_id).t         = SampleTime;
    data.(my_id).y         = Values(ind,:);
    data.(my_id).sigma     = StdDev(ind,:);
    data.(my_id).r         = r;
    data.(my_id).rep       = SampleReplicate;
  end
end
