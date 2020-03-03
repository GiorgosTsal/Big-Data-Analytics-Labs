function [mutM] = mutual(xV,maxtau,partitions)
% [mutM] = mutual(xV,maxtau,partitions)
% This function computes the mutual information of a given
% time series 'xV' for a number of lags 'tau'=0,...,'maxtau',
% where 'maxtau' is given. The probabilities are evaluated
% partitioning the domain and the number of partitions in one
% dimension is given by 'partitions'.
% The output is a two column matrix, the first column has the
% lags and the second is the mutual information.
% If 'partitions' is not specified the default is 16, and if
% 'maxtau' is also not specified the default is 'maxtau'=10.

n = length(xV);
if nargin==2
    partitions = round(sqrt(n/5));
elseif nargin==1
    partitions = round(sqrt(n/5));
    maxtau = 10;
end

h1V = zeros(partitions,1);  % for p(x(t+tau))
h2V = zeros(partitions,1);  % for p(x(t))
h12M = zeros(partitions,partitions);  % for p(x(t+tau),x(t))

% Normalise the data
xmin = min(xV);
[xmax,imax] = max(xV);
xV(imax) = xmax + (xmax-xmin)*10^(-10); % To avoid multiple exact maxima
yV = (xV-xmin)/(xmax-xmin);

arrayV = floor(yV*partitions)+1; % Array of partitions: 1,...,partitions
arrayV(imax) = partitions; % Set the maximum in the last partition

mutM = zeros(maxtau+1,2);
mutM(1:maxtau+1,1) = [0:maxtau]';
for tau=0:maxtau
    ntotal = n-tau;
    mutS = 0;
    for i=1:partitions
        for j=1:partitions
            h12M(i,j) = length(find(arrayV(tau+1:n)==i & arrayV(1:n-tau)==j));
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
    mutM(tau+1,2) = mutS;
end


