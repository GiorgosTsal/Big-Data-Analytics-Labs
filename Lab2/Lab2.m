%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file

%% Lab 2: Partial Cross-correlation on gene expression data 
alpha = 0.05;
K = 10;
npart = 3; % The number of variables in the conditional set Z when computing
           % the partial correlation r(X,Y|Z). 
           % If 0-> Z contains all the rest K-2 variables, 
           % else Z is a subset of 'npart' variables which are the most
           % correlated with the two variables X and Y.
rthresh = 0.2;
M = 100;
dattxt = '../data/MDma.xlsx';
maxwordlength = 15;
rng(1);

%% Read the data in the data matrix 'xM' with n rows (experiments) and m columns (genes)
[xM,txtM]=xlsread(dattxt);
geneC = txtM(2:end,1);
experC = txtM(1,2:end);
xM = xM';
[n,m]=size(xM);

%% Select subset of K genes (K<=m)
iV = randperm(m);
xM = xM(:,iV(1:K));
geneC = geneC(iV(1:K));

%% Compute the correlation matrix and the significance (p-values)
[ccM,p1M] = corrcoef(xM);
ccM(1:K+1:K*K) = 0; % assign diagonal values to zero
if K==2
    pccM = ccM;  % K=2: partial cross-correlation = cross-correlation
else
    pccM = NaN(K,K); 
    p1M = NaN(K,K); 
end
%% Compute partial cross-correlation for each pair of variables
if npart==0
    npart = K;
end
for i=1:K-1
    for j=i+1:K
        %% Compute the r(X_i,X_j|Z) but first define Z 
        xindV = setdiff([1:K],[i j]);
        if npart>=K-2
            indZV = xindV;
        else
            % Find the 'npart' most correlated variables to X_i and X_j
            cciV = ccM(i,xindV);
            ccjV = ccM(j,xindV);
            [occV,indoccV]=sort(abs(cciV)+abs(ccjV),'descend');
            indZV = xindV(indoccV(1:npart));
        end 
        % Find the regression residual of X_i to Z and the same for X_j
        [~,~,e1V] = regress(xM(:,i),[ones(n,1) xM(:,indZV)]);
        [~,~,e2V] = regress(xM(:,j),[ones(n,1) xM(:,indZV)]);
        tmpM = corrcoef(e1V,e2V);
        pccM(i,j)=tmpM(1,2);
        pccM(j,i)=pccM(i,j);
        %% Make a t-parametric test of significance for partial correlation 
        tstat = sqrt(n-2)*pccM(i,j)/sqrt(1-pccM(i,j).^2);
        p1M(i,j) = 2*(1-tcdf(abs(tstat),n-2));
        p1M(j,i)=p1M(i,j);
    end
end
tit1txt = sprintf('R_{XY|Z},#Z=%d',npart);
h1 = plotnetworktitle(abs(pccM),[],geneC,tit1txt,1,0,maxwordlength);

adj1M = p1M < alpha;
tit2txt = sprintf('Parametric p(R_{XY|Z}) < %1.2f',alpha);
h2 = plotnetworktitle(adj1M,[0 1],geneC,tit2txt,2,0,maxwordlength);

adjfdr1M = adjFDRmatrix(p1M,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f), parametric R_{XY|Z}',alpha);
h3 = plotnetworktitle(adjfdr1M,[0 1],geneC,tit3txt,3,0,maxwordlength);

adjthrM = abs(ccM) > rthresh;
tit4txt = sprintf('Adjacency R_{XY|Z} > %1.2f',rthresh);
h4 = plotnetworktitle(adjthrM,[0 1],geneC,tit4txt,4,0,maxwordlength);

