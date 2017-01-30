% Script replicate_regression_omics_selected
%
% Display results for a single protein; save graphics files.
% To run this script, a replicate regression for the entire omics data set must have been run before
%
% (see script replicate_regression_omics_analysis)
%
% Input variable: foptions_file
% Example:        foptions_file = [replicate_regression_PATH_NAME '/options_hahne_salt_stress_tagD.csv'];


% ----------------------------------------------------------------------------------------
% Options

foptions_default = join_struct(replicate_regression_omics_default_options, ...
                     struct('plot_data', 1,'plot_regression',1,'plot_replicates', 1,'plot_all',0,...
                          'graphics_format', 'png', 'graphics_scale','absolute','title_string',[],'delimiter_symbol','|',...
                          'x_label','Time','y_label','','convenience_name',[],'postprocess_normalise',0,...
                          'log_transformation','arithmetic','name_print_capital',0,'run_crossvalidation',...
                          0,'gp_subplot', []));

foptions = load_options_table(foptions_file);
foptions = join_struct(foptions_default,foptions);


% ----------------------------------------------------------------------------------------
% Load regression results

cd(foptions.data_dir);     load(foptions.data_file_matlab);
cd(foptions.result_dir);   load(foptions.result_file_matlab);
cd(foptions.graphics_dir);
ca;


% ----------------------------------------------------------------------------------------
% Initialise

replicates    = fieldnames(data_replicates);
n_rep         = length(replicates);
dum           = fieldnames(data_replicates.(replicates{1}));
element_names = getfield(data_reg.combined,(dum{1}));

if isfield(foptions,'element_id'),
  idstring = foptions.element_id; 
  if isnumeric(idstring), idstring = num2str(idstring); end
  id_list  = Strsplit(foptions.delimiter_symbol,idstring);
  if length(foptions.convenience_name),
    itt       = label_names(foptions.convenience_name,dum);
    ll        = label_names(id_list,data_reg.combined.(dum{itt}));
    name_list = data_reg.combined.(foptions.convenience_name)(ll);
  else, 
    name_list = id_list;
  end
else,
  namestring = foptions.element_name;
  name_list  = Strsplit(foptions.delimiter_symbol,namestring)';
  allnames = data_reg.combined.(foptions.convenience_name);
  %% the following two lines are there for safety reasons (in case names are missing)
  ind_empty = find(cellfun('isempty',allnames)); 
  allnames(ind_empty) = repmat({''},length(ind_empty),1);
  ll       = label_names(name_list,allnames);
  if find(ll==0), 
    name_list(find(ll==0))
    error('Unknown element label found');
  end
  id_list  = data_reg.combined.(dum{1})(ll);
end

if foptions.name_print_capital,
 for it = 1:length(name_list),
   name_list{it} = [upper(name_list{it}(1)) name_list{it}(2:end)];
 end
end


% ----------------------------------------------------------------------------------------
% Graphics options

colors  = {[0 0.3 1],[1 0.2 .2],[.9 .6 0],[0.3 1 0],[1 0.3 1]};

clear gp; 
gp.fignum          = 1; 
gp.subplot         = []; 
gp.replicate_names = fieldnames(data_reg.replicates); 
gp.fontsize        = 28; 
gp.x_label         = foptions.x_label; 
gp.y_label         = foptions.y_label; 
gp.no_legend             = 1;
gp.print_title     = 0;
gp.convenience_name = foptions.convenience_name; 
if strcmp(foptions.graphics_scale,'log2'),  
  if strcmp(foptions.data_scale,'absolute'), 
    gp.show_log2 = 1;
  else
    gp.logarithmic_data = 1;
  end; 
end
gp.image_format    = foptions.graphics_format;


% ----------------------------------------------------------------
% Graphics

ind_element = label_names(id_list,element_names);


%% Show data only
  
if foptions.plot_data,
  if length(data_pointwise_average),
    gp.data_lines      = 1;
    gp.flag_only_data  = 0;
    gp.flag_omit_replicates = 1;
    gp.title_string    = [foptions.title_string ' data'];
    biotable_interpolation_graphics_std(data_reg.combined, data_pointwise_average, data_reg.replicates, gp, ind_element, [ foptions.graphics_file '_data']);
  end
end

  
%% Show regression with replicate curves

if foptions.plot_replicates,
  gp.fignum          = 2; 
  gp.data_lines      = 0;
  gp.flag_only_data  = 0;
  gp.flag_omit_replicates = 0;
  gp.title_string    = [foptions.title_string ' replicates'];
  biotable_interpolation_graphics_std(data_reg.combined, data_reg.average, data_reg.replicates, gp, ind_element, [ foptions.graphics_file '_replicates']);
end


%% Show regression without replicate curves

if foptions.plot_regression,
  gp.fignum                = 3; 
  gp.flag_omit_replicates  = 1;
  gp.title_string    = [foptions.title_string ' regression'];
  biotable_interpolation_graphics_std(data_reg.combined, data_reg.average, data_reg.replicates, gp, ind_element, [ foptions.graphics_file '_regression']);

  if foptions.postprocess_normalise,
    if isfield(data_reg,'adjusted'),
      gp.fignum                = 4; 
      gp.data_lines            = 1;
      gp.flag_omit_replicates  = 1;
      gp.title_string    = [foptions.title_string ' normalised'];
      biotable_interpolation_graphics_std(data_reg.adjusted.combined, data_reg.adjusted.average, data_reg.adjusted.replicates, gp, ind_element, [ foptions.graphics_file '_normalised']);
    end
  end
  
end


%% Show all elements together (only data)

if foptions.plot_all,
   
  figure(5); clf; h = []; set(gca,'FontSize',12);
  cm = jet(length(ind_element));
  hold on;
  for it = 1:length(ind_element),

    %% Adjusted data
    %%h(it) = plot(data_reg.adjusted.combined.SampleTime',data_reg.adjusted.combined.DataMean(ind_element(it),:),'o','Color',cm(it,:));
    %%plot(data_reg.adjusted.average.SampleTime',data_reg.adjusted.average.DataMean(ind_element(it),:),'-','Color',cm(it,:));

    %% Non-adjusted data pointwise
    %%h(it) = plot(data_pointwise_average.SampleTime',data_pointwise_average.DataMean(ind_element(it),:),'o','Color',cm(it,:));
    %%plot(data_reg.average.SampleTime',data_reg.average.DataMean(ind_element(it),:),'-','Color',cm(it,:));

    %% Non-adjusted data
    h(it) = plot(data_reg.combined.SampleTime',data_reg.combined.DataMean(ind_element(it),:),'o','Color',cm(it,:));
    plot(data_reg.average.SampleTime',data_reg.average.DataMean(ind_element(it),:),'-','Color',cm(it,:));
  end
  legend(h,name_list);
  xlabel(foptions.x_label);
  ylabel('Protein level (a.u.)');
  title(foptions.title_string);

  cd( foptions.graphics_dir);
  switch gp.image_format,
    case 'png', print([foptions.graphics_file '.png'],'-f5','-dpng');
    case 'eps', print([foptions.graphics_file '.eps'],'-f5','-depsc');
  end

  if foptions.run_crossvalidation == 0,
    foptions.mark_outliers_percentage = 0;
  end
  
  if foptions.mark_outliers_percentage,
     %% deviation between crossvalidated and normal fit (replicate curves)
     deviation  = abs(data_reg.crossvalidation_replicate.DataMean - data_reg.fit.DataMean);
     %% deviation between crossvalidated fit and data point (replicate curves)
     qq         = quantile(deviation(:),1-foptions.mark_outliers_percentage);
     presumable_outliers = sparse(deviation>qq);  
   else
     presumable_outliers = data_reg.presumable_outliers;
  end
  
  figure(6); clf; 
  cd(foptions.graphics_dir);
  clear gp; 
  gp.fignum             = 6;
  gp.show_errorbars     = 1;
  gp.subplot            = [4,4]; 
  if length(foptions.gp_subplot),
    gp.subplot = foptions.gp_subplot;
  end
  gp.log_transformation = foptions.log_transformation; 
  gp.replicate_names    = fieldnames(data_reg.replicates); 
  gp.fontsize           = 10;
  gp.image_format       = foptions.graphics_format;
  gp.convenience_name   = foptions.convenience_name;
  gp.name_print_capital = foptions.name_print_capital;
  switch foptions.data_scale,
    case 'absolute', gp.logarithmic_data = 0; 
    otherwise,       gp.logarithmic_data = 1; 
  end 
  switch foptions.graphics_scale,
    case 'log2',   gp.show_log2 = 1;
    case 'absolute', gp.show_log2 = 0;
    otherwise,     gp.show_log2 = 0;
  end;
  gp.show_bands         = 1;
  gp.mark_data          = presumable_outliers;
  biotable_interpolation_graphics_std(data_reg.combined, data_reg.average, data_reg.replicates, gp, ind_element,[foptions.graphics_file '_panels']);
end


% Find the curve that is most distant from average curve: 

[dum, extremist] = max(sum([ data_reg.average.DataMean(ind_element,:) - repmat(mean(data_reg.average.DataMean(ind_element,:)), length(ind_element),1)].^2,2));

display(sprintf('Outermost item: %s',name_list{extremist}))
display(sprintf('Graphics directory: %s',foptions.graphics_dir));
