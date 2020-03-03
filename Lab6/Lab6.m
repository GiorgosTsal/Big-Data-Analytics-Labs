% Lab6: cross mutual information in greek stocks and networks
%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file

%% load data and set parameters
tau = 0;
alpha = 0.05;
K = 4;
nonparametric = 1; 
mithresh = 0.1;
M = 100;
maxtau = 20;

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

% Use log returns
xM = log(yM(2:n,:))-log(yM(1:n-1,:));

% For each pair compute the matrix of cross mutual information
miM = NaN*ones(K,K);
if tau==0
    % The matrix of cross mutual information is symmetric
    mutM = mutualxy(xM);
    count = 0;
    for ik=1:K-1
        for jk=ik+1:K
            count=count+1;
            miM(ik,jk)=mutM(count,3);
            miM(jk,ik)=mutM(count,3);
        end
    end
    if nonparametric
        % Randomized r values for the pair (y1,y2)
        p2M = NaN*ones(K,K);
        misurT = NaN*ones(M,K,K);
        misurT(1,:,:) = miM;
        for i=1:M
            zM = NaN*ones(n-1,K);
            for j=1:K
                randindx = ceil((n-maxtau)*rand);
                isurV = [randindx+1:n-1 1:randindx]';
                zM(:,j) = xM(isurV,j);
            end
            mutM = mutualxy(zM);
            count = 0;
            for ik=1:K-1
                for jk=ik+1:K
                    count=count+1;
                    misurT(i+1,ik,jk)=mutM(count,3);
                    misurT(i+1,jk,ik)=mutM(count,3);
                end
            end
        end
        for ik=1:K-1
            for jk=ik+1:K
                mixyV = squeeze(misurT(:,ik,jk));
            	p2M(ik,jk) = resampledpvalue(mixyV,1);
                p2M(jk,ik) = p2M(ik,jk);
            end
        end
    end
else
    % The correlation matrix is not symmetric
    for ik=1:K-1
        for jk=ik+1:K
            tmpM = mutualxy([xM(1:end-tau,ik) xM(1+tau:end,jk)]);
            miM(ik,jk) = tmpM(1,3);
            tmpM = mutualxy([xM(1:end-tau,jk) xM(1+tau:end,ik)]);
            miM(jk,ik) = tmpM(1,3);
        end
    end
    if nonparametric
        % Randomized r values for the pair (y1,y2)
        p2M = NaN*ones(K,K);
        misurT = NaN*ones(M,K,K);
        misurT(1,:,:) = miM;
        for i=1:M
            zM = NaN*ones(n-1,K);
            for j=1:K
                randindx = ceil((n-maxtau)*rand);
                isurV = [randindx+1:n-1 1:randindx]';
                zM(:,j) = xM(isurV,j);
            end
            % The correlation matrix is not symmetric
            for ik=1:K-1
                for jk=ik+1:K
                    tmpM = mutualxy([zM(1:end-tau,ik) zM(1+tau:end,jk)]);
                    misurT(i+1,ik,jk) = tmpM(1,3);
                    tmpM = mutualxy([zM(1:end-tau,jk) zM(1+tau:end,ik)]);
                    misurT(i+1,jk,ik) = tmpM(1,3);
                end
            end
        end
        for ik=1:K-1
            for jk=ik+1:K
                mixyV = squeeze(misurT(:,ik,jk));
            	p2M(ik,jk) = resampledpvalue(mixyV,1);
                mixyV = squeeze(misurT(:,jk,ik));
            	p2M(jk,ik) = resampledpvalue(mixyV,1);
            end
        end
    end
end
tit1txt = sprintf('I_{XY}(%d)',tau);
h1 = plotnetworktitle(miM,[],nameM,tit1txt,1);

adj2M = p2M < alpha;
tit2txt = sprintf('Adjacency p(R_{XY}(%d)) < %1.2f',tau,alpha);
h2 = plotnetworktitle(adj2M,[0 1],nameM,tit2txt,2);

adjfdr2M = adjFDRmatrix(p2M,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f), nonparametric, I_{XY}(%d)',alpha,tau);
h3 = plotnetworktitle(adjfdr2M,[],nameM,tit3txt,3);

mithreshM = miM > mithresh;
tit4txt = sprintf('Adjacency I_{XY}(%d) > %1.2f',tau,mithresh);
h4 = plotnetworktitle(mithreshM,[],nameM,tit4txt,4);

