function replicate_regression_single_display(t,y,sigma,t_true,x_true,result,options)

% replicate_regression_single_display(t,y,sigma,t_true,x_true,result,options)

options_default = struct('fignum',[1,2],'show_central',1,'fontsize',12,'show_sample_curves',0,'figure_title','');
eval(default('options','struct','t_true','[]','result','[]'));
options = join_struct(options_default,options);

figure(options.fignum(1)); clf; h = [];
set(gca,'FontSize',options.fontsize); 

if length(result),
  plot_range(result.t,result.x,[result.x-result.sigma; result.x+result.sigma],[],[0 0 1]); hold on;
end

h = []; 
legends = {};
if length(result.t),h = [h errorbar(t,y,sigma,'bo')];                 hold on; legends =[legends {'Data'}]; end
if length(t_true),  h = [h plot(t_true,x_true,'m--','Linewidth',2)];  legends = [legends {'True curve'}];   end

if length(h), legend(h,legends,'Location','NorthWest'); end
hold off;
title(options.figure_title);

% if options.show_sample_curves
%  figure(3); clf; set(gca,'FontSize',options.fontsize); hold on; h=[];
%  h(1) = plot(t_true,x_true,'b--','Linewidth',2); hold on;
%  h(2) = errorbar(t,y,sigma,'rs');  hold on;
%  h(5) = plot(result.t_int,result.x_sample_replicate(:,1),'r--');  hold on;
%  h(6) = plot(t_int,result.x_sample_replicate(:,2),'m--');  hold on;
%  h(7) = plot(t_int,result.x_sample_replicate(:,3),'g--');  hold on;
%  hold off;  
%  %axis([0 3 0 1.5])
%  legend(h,'True curve','Data 1','Data 2','Data 3','Sampled 1','Sampled 2','Sampled 3','Location','NorthWest');
% title(options.figure_title);
% end
