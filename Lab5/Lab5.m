%% Lab 5: cross correlations in greek stocks and networks
%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file

%% load data and set parameters
tau = 0;
alpha = 0.05;
K = 5;
rthresh = 0.2;
maxtau = 20;
dattxt = 'GRstocks';

yM = load('../data/stocks2003.dat');
[n,m]=size(yM);
% rng(1);

iV = randperm(m);
% iV = [1:m];
yM = yM(:,iV(1:K));

%% Read the names of the stocks
nameM = textread('../data/stock_names.dat','%s');
nameM = nameM(iV(1:K),:);

% If NaN replace them with interpolated values for each time series
for i=1:K
    i1V = find(isnan(yM(:,i)));
    if ~isempty(i1V)
        iokV = setdiff([1:n]',i1V);
        yM(i1V,i) = interp1(iokV,yM(iokV,i),i1V,'spline');
    end
end

%% Use log returns
xM = log(yM(2:n,:))-log(yM(1:n-1,:));

% For each pair compute the correlation matrix
ccM = NaN*ones(K,K);
p1M = zeros(K,K);
if tau==0
    % The correlation matrix is symmetric
    [ccM,p1M] = corrcoef(xM);
    p1M(1:K+1:K*K) = 0;
else
    % The correlation matrix is not symmetric
    for ik=1:K-1
        for jk=ik+1:K
            [tmpM,ptmpM] = corrcoef(xM(1:end-tau,ik),xM(1+tau:end,jk));
            ccM(ik,jk) = tmpM(1,2);
            p1M(ik,jk) = ptmpM(1,2);
            [tmpM,ptmpM] = corrcoef(xM(1:end-tau,jk),xM(1+tau:end,ik));
            ccM(jk,ik) = tmpM(1,2);
            p1M(jk,ik) = ptmpM(1,2);
        end
    end
end    
tit1txt = sprintf('R_{XY}(%d)',tau);
h1 = plotnetworktitle(ccM,[],nameM,tit1txt,1);

adj1M = p1M < alpha;
tit2txt = sprintf('Adjacency p(R_{XY}(%d)) < %1.2f',tau,alpha);
h2 = plotnetworktitle(adj1M,[0 1],nameM,tit2txt,2);

adjfdr1M = adjFDRmatrix(p1M,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f) R_{XY}(%d)',alpha,tau);
h3 = plotnetworktitle(adjfdr1M,[0 1],nameM,tit3txt,3);

rthreshM = abs(ccM) > rthresh;
tit4txt = sprintf('Adjacency R_{XY}(%d) > %1.2f',tau,rthresh);
h4 = plotnetworktitle(rthreshM,[0 1],nameM,tit4txt,4);

