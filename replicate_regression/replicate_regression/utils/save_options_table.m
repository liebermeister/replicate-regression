function my_table = save_options_table(options,tablefile)

% my_table = save_options_table(options, tablefile)
%
% Save matlab structure 'options' as a table file

fn = fieldnames(options);

for it = 1:length(fn),
  my_table{it,1} = fn{it};
  value = getfield(options,fn{it});
  if length(value) > 1,
    value = sprintf('%5.2f ',value); 
  end
  if isnumeric(value), value = num2str(value); end
  if isempty(value), value = ''; end
  my_table{it,2} = value;
end
 
if exist('tablefile','var'),
  mytable(my_table,0,tablefile);
end