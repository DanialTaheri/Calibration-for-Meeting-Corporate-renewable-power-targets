
% REC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
REC_steps=12;
load('NewJerseyREC');
PriceDates=nonzeros(SREC_price(:,1))+datenum([1900,01,01])-2;
Prices=nonzeros(SREC_price(:,2));
figure;
plot(PriceDates,Prices);
datetick();
title('NJ REC Prices');
xlabel('Date');
ylabel('Price ($)');
% Prices at t, X(t)
y=50;
X=Prices./y;
Y=X(2:end)-X(1:end-1);
Pt = X(2:end);
% Prices at t-1, X(t-1)
Pt_1 = X(1:end-1);
% Discretization for monthly prices
dt = 1/12;
count=0;
 
% PDF for discretized model
mrjpdf = @(Pt, sigmaSq, b, beta) ...
   betapdf(Pt,(2*b*beta)/sigmaSq, 2*(b)*(1-beta)/sigmaSq).* (ones(length(Pt),1)+f(sigmaSq, beta, b, Pt,Pt_1));
lb= [0.001 0.001 0.0001];
ub=[Inf Inf 0.9999];
t0 = [var(X)  0.03 0.17];
params= mle(Pt,'pdf',mrjpdf,'start',t0,'lowerbound',lb,'upperbound',ub,...
 'optimfun','fmincon');
acov = mlecov(params,Pt,'pdf',mrjpdf);
display(acov);
display(params);
% acov = mlecov(params,Pt,'logpdf',@betalogpdf);
SE = sqrt(diag(acov));
pVal=2*(tcdf(-abs(params'./SE),inf));
%print(pVal);
sigma=sqrt(params(1));
b= params(2);
beta=params(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Substitution in the paper
sigma_R=sqrt(params(1)*REC_steps);
Kappa= params(2)*REC_steps;
nu=beta*b*REC_steps;
% Simulate for about 3 years

 nPeriods = 40*12;
 nTrials = 1000;
 n1 = randn(nPeriods,nTrials);
 Simvalue = zeros(nPeriods, nTrials);
SimPrices=zeros(nPeriods, nTrials);
 Simvalue(1,:) = X(end);
 SimPriceDates=zeros(nPeriods,1);
 Y1=50;
 for i=2:nPeriods
     Simvalue(i,:) = Simvalue(i-1,:).*(1-b) + b*beta + ...
                 sigma*sqrt(Simvalue(i-1,:).*(1-Simvalue(i-1,:))).*n1(i,:);
             for j=1:nTrials
             if Simvalue(i,j)>=1
                 Simvalue(i,j)=0.99;
             end
             if Simvalue(i,j)<=0
                 Simvalue(i,j)=0.01;
                 count=count+1;
             end
             end
 end
 for i=1:nPeriods
 SimPrices(i,:)= Simvalue(i,:)*50;
 end
 for i=1:nPeriods
 SimPriceDates(i) = addtodate(PriceDates(end),i-1,'month');
 end
 figure;
 plot(SimPriceDates, SimPrices);
 title('New Jersey REC Prices simulation');
 xlabel('Date');
 ylabel('Price ($)');
 datetick();
 mean_SimPrices=zeros(nPeriods,1);
 for i=1:nPeriods
 mean_SimPrices(i,1)= sum(SimPrices(i,:))/nTrials;
 end
 s=zeros(length(Prices),nTrials);
s(1,:)= X(1);
n12=randn(length(Prices),nTrials);
Simvalue1=zeros(length(Prices),1);
Simvalue1(1,1)=X(1);
for i=2:length(Prices)
   s(i,:) = s(i-1,:)*(1-b) + b*beta + ...
                 sigma*sqrt(s(i-1,:).*(1-s(i-1,:))).*n12(i,:);
           for j=1:nTrials           
             if s(i,j)>=1
                 s(i,j)=0.99;
             end
             if s(i,j)<=0
                 s(i,j)=0.01;
             end
           end
           Simvalue1(i,1)=sum(s(i,:))/nTrials;
end
figure;
plot(PriceDates,Simvalue1);
hold on;
plot(PriceDates,X);
datetick();
hold off
figure;
plot(PriceDates,s(:,2));
hold on;
plot(PriceDates,X);
datetick();
hold off
title('New Jersey REC Prices simulation versus historical price');
xlabel('Date');
ylabel('Price ($)');
legend('Simulated price', 'Historical price');
figure;
qqplot(log(X(2:end)./X(1:end-1)), log(s(2:end,2)./s(1:end-1,2)));
title('QQ of Jacobi model');
xlabel('Quantile of REC return');
ylabel('Quantile of Jacobi diffusion');



