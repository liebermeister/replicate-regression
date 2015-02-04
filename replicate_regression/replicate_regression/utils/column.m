function x = column(x)

% transform a vector (row or column) into a column vector

if size(x,2)>1, x = x.'; end 