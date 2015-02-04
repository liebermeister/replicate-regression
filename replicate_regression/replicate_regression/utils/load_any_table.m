function result = load_any_table(filename,delimiter)

% list = load_any_table(filename)
%
% loads a tab-delimited file and puts it into a cell array
% filename: name of tab-delimited file containing strings and numbers

if ~exist('delimiter','var'), delimiter = sprintf('\t'); end

fid           = fopen(filename);
try
  column_titles = fgets(fid);
catch
  error(['File ' filename ' not found.']);
end
fclose(fid);

tab_pos = [0,unique([strfind(column_titles,delimiter)])];

textscanstring = repmat('%s',1,length(tab_pos));

fid = fopen(filename);
A   = textscan(fid,textscanstring,'delimiter',delimiter);
fclose(fid);

result = {};
for i =1:length(A),
  for k = 1:length(A{i}),
  result{k,i} = A{i}{k};
  end
end

result = result(:,find(sum(cellfun('length',result))>0));
