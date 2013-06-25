function my_table = save_options_table(options,tablefile)

fn = fieldnames(options);

for it = 1:length(fn),
  my_table{it,1} = fn{it};
  value = getfield(options,fn{it});
   if isnumeric(value), value = num2str(value); end
  my_table{it,2} = value;
end

if exist('tablefile','var'),
  table(my_table,0,tablefile);
end