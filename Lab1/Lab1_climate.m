%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file
%% Lab 1: Cross-correlation on gene expression data 
alpha = 0.05;
K = 10;
nonparametric = 1; 
rthresh = 0.2;
M = 100;
dattxt = '../data/MDma.xlsx';
dattxt2 = '../data/jena_climate_2009_2016.csv';
maxwordlength = 15;
rng(1);

%% Read the data in the data matrix 'xM' with n rows (experiments) and m columns (genes)
[xM,txtM]=xlsread(dattxt);
geneC = txtM(2:end,1);
experC = txtM(1,2:end);
xM = xM';
[n,m]=size(xM);


%% my way
dataM = readtable(dattxt2);
%% Select subset of K genes (K<=m)
iV = randperm(m);
xM = xM(:,iV(1:K));
geneC = geneC(iV(1:K));
%% Convert first column(datetime) to number-not working
%dataM.DateTime = datenum(dataM.DateTime, 'dd-mm-yyyy HH:MM:SS');
%% Plot all variables together
figure(20);
stackedplot(dataM);
%% Compute the correlation matrix and the significance (p-values)
[ccM,p1M] = corrcoef(xM);
p1M(1:K+1:K*K) = 0; % assign diagonal values to zero
if nonparametric
    % Randomized r values for the pair (x1,x2)
    p2M = zeros(K,K);
    ccsurT = NaN*ones(M,K,K);
    ccsurT(1,:,:) = ccM;
    for i=1:M
        zM = NaN*ones(n,K);
        for j=1:K
            zM(:,j) = xM(randperm(n),j);
        end
        ccsurT(i+1,:,:) = corrcoef(zM);
    end
    for ik=1:K-1
        for jk=ik+1:K
            rxyV = squeeze(ccsurT(:,ik,jk));
            p2M(ik,jk) = resampledpvalue(rxyV,1);
            p2M(jk,ik) = p2M(ik,jk);
        end
    end
end
tit1txt = sprintf('R_{XY}');
h1 = plotnetworktitle(abs(ccM),[],geneC,tit1txt,1,0,maxwordlength);

adj1M = p1M < alpha;
tit2txt = sprintf('Parametric p(R_{XY}) < %1.2f',alpha);
h2 = plotnetworktitle(adj1M,[0 1],geneC,tit2txt,2,0,maxwordlength);

adjfdr1M = adjFDRmatrix(p1M,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f), parametric R_{XY}',alpha);
h3 = plotnetworktitle(adjfdr1M,[0 1],geneC,tit3txt,3,0,maxwordlength);

adjthrM = abs(ccM) > rthresh;
tit4txt = sprintf('Adjacency R_{XY} > %1.2f',rthresh);
h4 = plotnetworktitle(adjthrM,[0 1],geneC,tit4txt,4,0,maxwordlength);

if nonparametric
    adj2M = p2M < alpha;
    tit5txt = sprintf('Randomization p(R_{XY}) < %1.2f',alpha);
    h5 = plotnetworktitle(adj2M,[0 1],geneC,tit5txt,5,0,maxwordlength);
    adjfdr2M = adjFDRmatrix(p2M,alpha,2);
    tit6txt = sprintf('FDR (a=%1.3f), randomization R_{XY}',alpha);
    h6 = plotnetworktitle(adjfdr2M,[0 1],geneC,tit6txt,6,0,maxwordlength);
end
