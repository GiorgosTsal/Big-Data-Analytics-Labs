function prcV = percentile(xM,percent)
% prcV = percentile(xM,percent)
% PERCENTILE finds the given percentile of a given list (vector). 
% If the input is a matrix it does this for each column of the matrix.
% INPUT 
% - xM      : the vector for which the percentile is computed. If 'xM' is a
%             matrix it is repeated for each column.
% - percent : the percentage for which the percentile is computed.
% OUTPUT 
% - prcV    : the percentile if input is a vector or the row vector of
%             percentiles if input is a matrix.

[n,m]=size(xM);
iprc = round(percent*n/100);
prcV = NaN(1,m);
for i=1:m
    xV = sort(xM(:,i));
    prcV(i) = xV(iprc);
end
