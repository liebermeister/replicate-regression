function h = replicate_regression_display(t,y,sigma,l,t_true,x_true,result,options)

% h = replicate_regression_display(t,y,sigma,l,t_true,x_true,result,options)
%
% Graphical output for results of replicate time series regression

options_default = struct('fignum',[1,2],'show_central',1,'fontsize',12,'show_sample_curves',0,'figure_title','','show_crossvalidation',1,'show_bands',1,'show_errorbars',1);

eval(default('options','struct','t_true','[]','result','[]','linewidth','1'));

options = join_struct(options_default,options);

colors = rr_colors;


% ---------------------------------
% Display results

t1     = t(l==1);
t2     = t(l==2);
t3     = t(l==3);
y1     = y(l==1);
y2     = y(l==2);
y3     = y(l==3);
sigma1 = sigma(l==1);
sigma2 = sigma(l==2);
sigma3 = sigma(l==3);

if length(result),
  t_reg            = result.t;
  if options.show_central,
    x_central        = result.x_central;
    sigma_central    = result.sigma_central;
  end
  x_reg            = result.x_average;
  sigma_reg        = result.sigma_average;
  x_reg_replicate     = result.x_replicate;
  sigma_reg_replicate = result.sigma_replicate;
else
  sigma_reg_replicate = [];
end


figure(options.fignum(1)); clf; h = [];
set(gca,'FontSize',options.fontsize); 
h = []; legends = {}; 


if options.show_central,
  if length(result),
    if ~options.show_bands,
      h = plot_range(t_reg,x_reg,[x_reg-sigma_reg; x_reg+sigma_reg],[],[0 0 0]); hold on;
    end
    legends = [legends {'Regression'}]; 
  end
end

if length(sigma_reg_replicate),
if options.show_bands,
  plot_range(t_reg,x_reg_replicate{1},[x_reg_replicate{1}-sigma_reg_replicate{1}; x_reg_replicate{1}+sigma_reg_replicate{1}],[],colors{1},[],[],[],'-',2); hold on;
  plot_range(t_reg,x_reg_replicate{2},[x_reg_replicate{2}-sigma_reg_replicate{2}; x_reg_replicate{2}+sigma_reg_replicate{2}],[],colors{2},[],[],[],'-',2); hold on;
  plot_range(t_reg,x_reg_replicate{3},[x_reg_replicate{3}-sigma_reg_replicate{3}; x_reg_replicate{3}+sigma_reg_replicate{3}],[],colors{3},[],[],[],'-',2);hold on;
else,
  plot(t_reg,x_reg_replicate{1},'-','Color',colors{1},'Linewidth',2); hold on;
  plot(t_reg,x_reg_replicate{2},'-','Color',colors{2},'Linewidth',2); hold on;
  plot(t_reg,x_reg_replicate{3},'-','Color',colors{3},'Linewidth',2); hold on;
end
end

if options.show_central,
  if length(result),
    if options.show_bands,
      h = [h plot(t_reg,x_reg,'k-','Linewidth',2)]; 
      plot(t_reg,x_reg + sigma_reg,'k--');  
      plot(t_reg,x_reg - sigma_reg,'k--');  
      legends = [legends {'Regression'}]; 
    end
  end
end


hold on;

if options.show_errorbars,
  if length(t1), h = [h errorbar(t1,y1,sigma1,'.','MarkerSize',20,'Color',colors{1})];  legends =[legends {'Data 1'}]; end
  if length(t2), h = [h errorbar(t2,y2,sigma2,'.','MarkerSize',20,'Color',colors{2})];  legends =[legends {'Data 2'}]; end
  if length(t3), h = [h errorbar(t3,y3,sigma3,'.','MarkerSize',20,'Color',colors{3})];  legends =[legends {'Data 3'}]; end
else,
  if length(t1), h = [h plot(t1,y1,'.','MarkerSize',20,'Color',colors{1})];  legends =[legends {'Data 1'}]; end
  if length(t2), h = [h plot(t2,y2,'.','MarkerSize',20,'Color',colors{2})];  legends =[legends {'Data 2'}]; end
  if length(t3), h = [h plot(t3,y3,'.','MarkerSize',20,'Color',colors{3})];  legends =[legends {'Data 3'}]; end
end


if length(t_true), h = [h plot(t_true,x_true,'m--','Linewidth',2)];  legends = [legends {'True curve'}];   end

if options.show_crossvalidation,
  if isfield(result,'x_crossvalidation');
    xc1     = result.x_crossvalidation(l==1);
    xc2     = result.x_crossvalidation(l==2);
    xc3     = result.x_crossvalidation(l==3);
    h = [h plot(t1,xc1,'*','Color',colors{1})];  legends = [legends {'Crossvalidation 1'}]; 
    h = [h plot(t2,xc2,'*','Color',colors{2})];  legends = [legends {'Crossvalidation 2'}]; 
    h = [h plot(t3,xc3,'*','Color',colors{3})];  legends = [legends {'Crossvalidation 3'}]; 
  end
end

if length(h), legend(h,legends,'Location','Best'); end
hold off; title(options.figure_title);

axis tight

if options.show_sample_curves
 figure(3); clf; set(gca,'FontSize',options.fontsize);  h=[];
 h(1) = plot(t_true,x_true,'b--','Linewidth',2); 
 h(2) = errorbar(t1,y1,sigma1,'rs');  
 h(3) = errorbar(t2,y2,sigma2,'ms');  
 h(4) = errorbar(t3,y3,sigma3,'gs');  
 h(5) = plot(t_reg,result.x_sample_replicate(:,1),'r--');  
 h(6) = plot(t_reg,result.x_sample_replicate(:,2),'m--');  
 h(7) = plot(t_reg,result.x_sample_replicate(:,3),'g--');  
 %axis([0 3 0 1.5])
 hold off;  
 legend(h,'True curve','Data 1','Data 2','Data 3','Sampled 1','Sampled 2','Sampled 3','Location','NorthWest');
title(options.figure_title);
end

