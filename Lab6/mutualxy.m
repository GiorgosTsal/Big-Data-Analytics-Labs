function [mutM] = mutualxy(xM,partitions)
% [mutM] = mutualxy(xM,partitions)
% This function computes the mutual information between variables
% given as columns of the input matrix 'xM'. If size(xM,1) > 2, all  
% combinations of the variables (columns) are considered. 
% The probabilities are evaluated partitioning the domain and 
% the number of partitions in one dimension is given by 'partitions'.
% The output is a matrix of  for the mutual information. 
% If 'partitions' is not specified the default is 16. 
% INPUT:
% xM        : sample matrix of all variables (one per column)
% paritions : the number of bins to split the domain of each variable
% OUTPUT:
% mutM      : 3 column matrix:
%             first column the index number of the first variable,
%             second column the index number of the second variable,
%             third column the mutual information of the two variables

[n,m] = size(xM);
if nargin==1
    partitions =  fix(sqrt(n/5));
end
ncomb = factorial(m)/(2*factorial(m-2));

% Normalise each column of the matrix 
xminV = min(xM);
[xmaxV,imaxV] = max(xM);
for j=1:m
    xM(imaxV(j),j) = xmaxV(j) + (xmaxV(j)-xminV(j))*10^(-10); % To avoid multiple exact maxima
end
yM = (xM-ones(n,1)*xminV)./(ones(n,1)*(xmaxV-xminV));
arrayM = floor(yM*partitions)+1; % Array of partitions: 1,...,partitions
for j=1:m
    arrayM(imaxV(j),j) = partitions; % Set the maximum in the last partition
end

h1V = zeros(partitions,1);  % for p(x)
h2V = zeros(partitions,1);  % for p(y)
h12M = zeros(partitions,partitions);  % for p(x,y)
mutM = NaN*ones(ncomb,3);
count = 0;
for i1=1:m-1
    for i2=i1+1:m
        count = count+1;
        mutS = 0;
        for i=1:partitions
            for j=1:partitions
                h12M(i,j) = length(find(arrayM(:,i1)==i & arrayM(:,i2)==j));
            end
        end
        for i=1:partitions
            h1V(i) = sum(h12M(i,:));
            h2V(i) = sum(h12M(:,i));
        end
        for i=1:partitions
            for j=1:partitions
                if h12M(i,j) > 0
                    mutS=mutS+(h12M(i,j)/n)*log(h12M(i,j)*n/(h1V(i)*h2V(j)));
                end
            end
        end
        mutM(count,1) = i1;
        mutM(count,2) = i2;
        mutM(count,3) = mutS;
    end
end
