function biotable_interpolation_graphics_std(combined,averaged,replicate_averaged,p,ind_show,outfile);

eval(default('outfile','[]','replicate_names','[]'));

p_default = struct('fignum',1, ...
		   'show_log2',0,...
		   'logarithmic_data',0,...
                   'log_transformation', 'arithmetic', ...
		   'subplot',[],...
                   'show_errorbars',1,...
		   'fontsize',12,...
                   'show_bands', 1, ...
		   'flag_only_data',0,...
		   'flag_omit_replicates',0,...
		   'linewidth',1,...
		   'markerwidth',1,...
		   'image_format','eps',...
		   'show_labels',1, ...
                   'data_lines', 0, ...
                   'title_string', [], ...
                   'print_title', 1, ...
                   'name_print_capital', 0, ...
                   'x_label', '',...
                   'y_label', '',...
                   'no_legend', 0,...
                   'convenience_name', [], ...
                   'mark_data', []);

p  = join_struct(p_default,p);

pname = [];

%% if replicates are stored in fields -> reformat them as a list

if isstruct(replicate_averaged),
  dum = replicate_averaged; 
  replicate_averaged = {};
  fn = fieldnames(dum);
  for it = 1:length(fn)
    replicate_averaged{it} = getfield(dum,fn{it});
  end
end

if p.show_log2,
  if ~p.logarithmic_data,
    [combined.DataMean, combined.DataStd] = lognormal_normal2log(combined.DataMean, combined.DataStd,p.log_transformation);
    combined.DataMean = combined.DataMean/log(2);
    combined.DataStd  = combined.DataStd/log(2);
    [averaged.DataMean, averaged.DataStd] = lognormal_normal2log(averaged.DataMean, averaged.DataStd,p.log_transformation);
    averaged.DataMean = averaged.DataMean/log(2);
    averaged.DataStd  = averaged.DataStd/log(2);
    for it = 1:length(replicate_averaged),
      [replicate_averaged{it}.DataMean, replicate_averaged{it}.DataStd] = ...
          lognormal_normal2log(replicate_averaged{it}.DataMean, replicate_averaged{it}.DataStd,p.log_transformation);
      replicate_averaged{it}.DataMean = replicate_averaged{it}.DataMean/log(2);
      replicate_averaged{it}.DataStd  = replicate_averaged{it}.DataStd/log(2);
  end
  end
end

n_rep = length(replicate_averaged);

if ~isfield(combined,'ProteinName'), 
  if isfield(combined,'GeneName'), 
    combined.ProteinName = combined.GeneName; 
  elseif isfield(combined,'MetaboliteName'), 
    combined.ProteinName = combined.MetaboliteName; 
  else
    dum = fieldnames(combined); 
    combined.ProteinName = getfield(combined,dum{1});
  end
end

eval(default('ind','1:length(combined.ProteinName)'));
eval(default('ind_show','ind'));

colors = rr_colors;

%if outfile, 
%  sfigure(p.fignum); %clf;
%end

for itt = 1:length(ind_show),

  it = ind_show(itt);
  
  if length(p.subplot)==0,
    sfigure(itt); clf;
  else,
    sfigure(p.fignum); 
    subplot(p.subplot(1),p.subplot(2),mod(itt-1,prod(p.subplot))+1); 
    p.replicate_names = [];
  end

  x_mean  = averaged.DataMean(it,:);
  x_std   = averaged.DataStd(it,:);
  x_upper = x_mean + x_std;
  x_lower = x_mean - x_std;
  h       = [];
  
  hold on;

  if sum(isfinite(x_mean)),

    %% average curve
    
    if ~p.flag_only_data,
      if p.show_bands ==0,
        h(1) = plot_range(column(averaged.SampleTime),column(x_mean),column(x_std),[],[0 0 0]);        hold on
      end
    end

    %% interpolated curves

    if ~[p.flag_only_data+p.flag_omit_replicates],
      for it_r = 1:n_rep,
        ind     = find(it_r==cell2mat(combined.SampleName));
        present = sum(isfinite(combined.DataMean(it,ind)));
        if present,
          if p.show_bands,
            plot_range(column(replicate_averaged{it_r}.SampleTime),column(replicate_averaged{it_r}.DataMean(it,:)),...
                       column(replicate_averaged{it_r}.DataStd(it,:)),[],colors{it_r});  hold on; 
          else
            plot(column(replicate_averaged{it_r}.SampleTime),column(replicate_averaged{it_r}.DataMean(it,:)),'-','Color',colors{it_r});  hold on; 
          end
        end
      end
    end
 
    %% average curve
    
    if ~p.flag_only_data,
      if p.flag_omit_replicates,
        h(1) = plot_range(column(averaged.SampleTime),column(x_mean),column(x_std),[],[0 0 0]);        hold on
      else,      
        h(1) = plot(column(averaged.SampleTime),column(x_mean),'-','Color',[0 0 0],'Linewidth',p.linewidth);  hold on;
        plot(column(averaged.SampleTime),column(x_upper),'--','Color', [0 0 0],'Linewidth',p.linewidth); hold on;
        plot(column(averaged.SampleTime),column(x_lower),'--','Color', [0 0 0],'Linewidth',p.linewidth); hold on
      end
    end

    %% data points

    for it_r = 1:n_rep, 
      ind     = find([it_r==cell2mat(combined.SampleName)] .* isfinite(combined.DataMean(it,:)'));
      present = sum(isfinite(combined.DataMean(it,ind)));
      if present,
        if p.show_errorbars,
          h(1+it_r) = errorbar(combined.SampleTime(ind),column(combined.DataMean(it,ind)),...
                               column(combined.DataStd(it,ind)),'.','Color',colors{it_r},...
                               'Linewidth',p.markerwidth,'Markersize',15*p.markerwidth); 
        else
          h(1+it_r) = plot(combined.SampleTime(ind),column(combined.DataMean(it,ind)),...
                           '.','Color',colors{it_r},'Markersize',15*p.markerwidth); 
        end
        hold on
        if p.data_lines,
          plot(combined.SampleTime(ind),column(combined.DataMean(it,ind)),...
               '-','Color',colors{it_r},'Linewidth',p.markerwidth);           hold on
        end
      else,
        h(1+it_r) = plot(0,nan,'.','Color',colors{it_r}); hold on
      end
      hold on
    end

    if length(p.mark_data),
      for it_r = 1:n_rep, 
        ind     = find([it_r==cell2mat(combined.SampleName)] .* isfinite(combined.DataMean(it,:)'));
        my_ind = ind(find(p.mark_data(it,ind)));
        plot(combined.SampleTime(my_ind),column(combined.DataMean(it,my_ind)),...
             'ok','Markersize',10*p.markerwidth);       
      end
    end
    
    %% labels etc
    %hold off; 
    set(gca,'FontSize',p.fontsize,'LineWidth',p.linewidth);
    
    if p.show_log2 + p.logarithmic_data,
      %%      set(gca,'YTickLabel',2.^cell_string2num(cellstr(get(gca,'YTickLabel'))));
       axis tight;
    else,
      axis_tight_positive; 
    end
    
    if p.show_labels,

      if length(p.convenience_name), 
          ts = combined.(p.convenience_name){it}; 
        else,
          ts = combined.ProteinName{it};
        end
        
        if p.name_print_capital,
          ts = [upper(ts(1)) ts(2:end)];
        end

        if p.print_title,
          tss = ts;
          if length(p.title_string), tss = [p.title_string ' ' ts]; end
          title(tss,'FontSize',p.fontsize);
        end
      
      if ~p.no_legend,
        if length(p.replicate_names),
          if ~p.flag_only_data,
            legend(h,[{'Mean'}; column(p.replicate_names)],'Location','Best'); 
          else,
            legend(h(2:end),[column(p.replicate_names)],'Location','Best'); 
          end 
        end
      end
    end
    
  end
 
  xlabel(p.x_label);   ylabel(p.y_label);
    
  %% save graphics
  
  if isempty(pname),
    pname = strrep(strrep(ts,' ','_'),'/','_');
  end
  
  if outfile,
    
    if length(p.subplot)==0,
      switch p.image_format,
        case 'eps',  print([ outfile '_' pname '.eps'],['-f' num2str(itt)],'-depsc');
          display(sprintf('Saving graphics file %s',[ outfile '_' pname '.eps']));
        case 'png',  print([ outfile '_' pname '.png'],['-f' num2str(itt)],'-dpng');
          display(sprintf('Saving graphics file %s',[ outfile '_' pname '.png']));
      end
      
    elseif [mod(itt,prod(p.subplot))==0] + [itt == length(ind_show)],
      
      switch p.image_format,
        case 'eps',  print([ outfile '_' pname '.eps'],['-f' num2str(p.fignum)],'-depsc');
          display(sprintf('Saving graphics file %s',[ outfile '_' pname '.eps']));
        case 'png',  print([ outfile '_' pname '.png'],['-f' num2str(p.fignum)],'-dpng');
          display(sprintf('Saving graphics file %s',[ outfile '_' pname '.png']));
      end
      clf
      pname = [];

    end
    
  end
  
end
