% Lab 4: Explore nonlinear autocorrelations and cross correlations in greek stocks 

%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file
%% Read data and set parameters
maxtau = 20;
maxtau2 = 10;
p = 1;
alpha = 0.05;
M = 100;

yM = load('../data/stocks2003.dat');
[n,m]=size(yM);
% rng(1);

% Read the names of the stocks
nameM = textread('../data/stock_names.dat','%s');

%% First index
indx1 = ceil(m*rand);
y1V = yM(:,indx1);
name1 = cell2mat(nameM(indx1,:));
% If NaN replace them with interpolated values
i1V = find(isnan(y1V));
if ~isempty(i1V)
    iokV = setdiff([1:n]',i1V);
    y1V(i1V) = interp1(iokV,y1V(iokV),i1V,'spline');
end

figure(1)
clf
plot(y1V,'.-')
xlabel('day t')
ylabel('y(t)')
title(sprintf('stock %s',name1))

miy1M = mutual(y1V,maxtau);
figure(2)
clf
plot(miy1M(:,1),miy1M(:,2),'.-')
xlabel('lag \tau')
ylabel('I_Y(\tau)')
title(sprintf('Delayed mutual information of stock %s',name1))

%% Second index
indx2 = ceil(m*rand);
y2V = yM(:,indx2);
name2 = cell2mat(nameM(indx2,:));
% If NaN replace them with interpolated values
i2V = find(isnan(y2V));
if ~isempty(i2V)
    iokV = setdiff([1:n]',i2V);
    y2V(i2V) = interp1(iokV,y2V(iokV),i2V,'spline');
end

figure(3)
clf
plot(y2V,'.-')
xlabel('day t')
ylabel('y(t)')
title(sprintf('stock %s',name2))

miy2M = mutual(y2V,maxtau);
figure(4)
clf
plot(miy2M(:,1),miy2M(:,2),'.-')
xlabel('lag \tau')
ylabel('I_Y(\tau)')
title(sprintf('Delayed mutual information of stock %s',name2))

mutyM = mutualxyts(y1V,y2V,maxtau2);
figure(5)
clf
plot(mutyM(:,1),mutyM(:,2),'.-')
xlabel('lag \tau')
ylabel('I_{XY}(\tau)')
title(sprintf('Cross mutual information of stock %s and %s',name1,name2))

%% Take log returns for the first stock index time series
ey1V = log(y1V(2:n))-log(y1V(1:n-1));
figure(6)
clf
plot(ey1V,'.-')
xlabel('day t')
ylabel('log returns of y(t)')
title(sprintf('stock %s, log returns',name1))

miey1M = mutual(ey1V,maxtau);
figure(7)
clf
plot(miey1M(:,1),miey1M(:,2),'.-')
xlabel('lag \tau')
ylabel('I_Y(\tau)')
title(sprintf('Delayed mutual information of stock returns %s',name1))

%% Randomized MI values for y1
misurM = NaN*ones(M,maxtau+1); 
for i=1:M
    iV = randperm(n-1);
    zV = ey1V(iV);
    tmpM = mutual(zV,maxtau);
    misurM(i,:) = tmpM(:,2)';
end
misurM = sort(misurM);
misur1V = percentile(misurM,(1-alpha)*100);

figure(7)
hold on
plot([0:maxtau],misur1V,'--c')

%% Take log returns for the second stock index time series
ey2V = log(y2V(2:n))-log(y2V(1:n-1));
figure(8)
clf
plot(ey2V,'.-')
xlabel('day t')
ylabel('log returns of y(t)')
title(sprintf('stock %s, log returns',name2))

miey2M = mutual(ey2V,maxtau);
figure(9)
clf
plot(miey2M(:,1),miey2M(:,2),'.-')
xlabel('lag \tau')
ylabel('I_Y(\tau)')
title(sprintf('Delayed mutual information of stock returns %s',name2))

%% Randomized MI values for y2
misurM = NaN*ones(M,maxtau+1); 
for i=1:M
    iV = randperm(n-1);
    zV = ey2V(iV);
    tmpM = mutual(zV,maxtau);
    misurM(i,:) = tmpM(:,2)';
end
misurM = sort(misurM);
misur2V = percentile(misurM,(1-alpha)*100);

figure(9)
hold on
plot([0:maxtau],misur2V,'--c')

mutxM = mutualxyts(ey1V,ey2V,maxtau2);
figure(10)
clf
plot(mutxM(:,1),mutxM(:,2),'.-')
xlabel('lag \tau')
ylabel('I_{XY}(\tau)')
title(sprintf('Cross mutual information of stock returns %s and %s',name1,name2))

%% Randomized MI values for the pair (y1,y2)
misurM = NaN*ones(M,2*maxtau2+1); 
for i=1:M
    randindx = ceil((n-maxtau)*rand);
    iV = [randindx+1:n-1 1:randindx]';
    zV = ey2V(iV);
    tmpM = mutualxyts(ey1V,zV,maxtau2);
    misurM(i,:) = tmpM(:,2)';
end
misurM = sort(misurM);
misurxyV = percentile(misurM,(1-alpha)*100);

figure(10)
hold on
plot([-maxtau2:maxtau2],misurxyV,'--c')

