function [mutM] = mutualxyts(xV,yV,maxtau,partitions)
% [mutM] = mutualxyts(xV,yV,maxtau,partitions)
% This function computes the cross delayed mutual information for two given
% time series 'xV' and 'yV' and for delays from -maxtau to maxtau.
% The probabilities are evaluated partitioning the domain of each variable
% (time series), and the number of partitions in each dimension is given by
% 'partitions'.
% The output is a two column matrix, the first column has the lags and the
% second is the mutual information.
% If 'partitions' is not specified the default is 16, and if
% 'maxtau' is also not specified the default is 'maxtau'=10.
% INPUT:
% xV        : a column vector of the first time series
% yV        : a column vector of the second time series
% paritions : the number of bins to split the domain of each variable
% OUTPUT:
% mutM      : 2 column matrix:
%             first column the lags (from -maxtau to maxtau)
%             second column the mutual information of the two time series

n = length(xV);
if nargin==3
    partitions = round(sqrt(n/5));
elseif nargin==2
    partitions = round(sqrt(n/5));
    maxtau = 10;
end

h1V = zeros(partitions,1);  % for p(x(t+tau))
h2V = zeros(partitions,1);  % for p(x(t))
h12M = zeros(partitions,partitions);  % for p(x(t+tau),x(t))

% Normalise both time series
xmin = min(xV);
[xmax,imax] = max(xV);
xV(imax) = xmax + (xmax-xmin)*10^(-10); % To avoid multiple exact maxima
xV = (xV-xmin)/(xmax-xmin);
ymin = min(yV);
[ymax,imax] = max(yV);
yV(imax) = ymax + (ymax-ymin)*10^(-10); % To avoid multiple exact maxima
yV = (yV-ymin)/(ymax-ymin);

arrayxV = floor(xV*partitions)+1; % Array of partitions: 1,...,partitions
arrayxV(imax) = partitions; % Set the maximum in the last partition
arrayyV = floor(yV*partitions)+1; % Array of partitions: 1,...,partitions
arrayyV(imax) = partitions; % Set the maximum in the last partition

mutM = zeros(2*maxtau+1,2);
mutM(1:2*maxtau+1,1) = [-maxtau:maxtau]';
for tau=-maxtau:maxtau
    ntotal = n-abs(tau);
    mutS = 0;
    for i=1:partitions
        for j=1:partitions
            if tau<0
                h12M(i,j) = length(find(arrayxV(-tau+1:n)==i & arrayyV(1:n+tau)==j));
            else
                h12M(i,j) = length(find(arrayxV(1:n-tau)==i & arrayyV(tau+1:n)==j));
            end                
        end
    end
    for i=1:partitions
        h1V(i) = sum(h12M(i,:));
        h2V(i) = sum(h12M(:,i));
    end
    for i=1:partitions
        for j=1:partitions
            if h12M(i,j) > 0
                mutS=mutS+(h12M(i,j)/ntotal)*log(h12M(i,j)*ntotal/(h1V(i)*h2V(j)));
            end
        end
    end
    mutM(tau+maxtau+1,2) = mutS;
end
