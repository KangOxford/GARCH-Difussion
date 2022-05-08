%% Option Pricing-GARCH





%% 1. Monte Carlo Simulation--Option Pricing (Model 1)   
% For Figure 1.
clear,clc
unit_t=1/365; r_1=0.039; 
mu=0.01; r=exp(r_1*unit_t); rho=-(mu-r_1*unit_t);

% ConVar(t)=Constant+GARCH*ConVar(t-1)+ARCH*Invt(t-1)^2; GARCH(1,1) Model;
alpha_0=0.01; alpha_1=0.5; beta=0.15;
Mdl = garch('Constant',alpha_0,'GARCH',beta,'ARCH',alpha_1);
[ConVar,Invt] = simulate(Mdl,350,'NumPaths',1);
z=Invt./sqrt(ConVar);

R=mu-0.5*ConVar+Invt;
S_0=1.5;
S=zeros(length(R),1);
S(1)=S_0*exp(R(1));
for t=2:length(R)
    S(t)=S(t-1)*exp(R(t));    
end

% Pricing below
t=100; T=102;
MC=100000;
u=randn(MC,T-t);

x1=zeros(MC,T-t);
x2=zeros(MC,T-t);
ConVF1=zeros(MC,T-t);
ConVF2=zeros(MC,T-t);
for j=1:T-t
    for i=1:MC
    ConVF1(i,j)=ConVolFun(x1(i,1:j),j,alpha_0,alpha_1,beta,ConVar(t),z(t));
    x1(i,j)=u(i,j)+rho/ConVF1(i,j)+ConVF1(i,j);
    
    ConVF2(i,j)=ConVolFun(x2(i,1:j),j,alpha_0,alpha_1,beta,ConVar(t),z(t));
    x2(i,j)=u(i,j)+rho/ConVF2(i,j);
    end
end

K=S(t)*1.25;
SL1=sum(ConVF1.*u,2);
SL2=sum(ConVF2.*u,2);
SR1=log(S(t)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(t)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC
Phi_D2=count2/MC
C=S(t)*Phi_D1-K*r^(t-T)*Phi_D2

% What is the shape of integral domains D_1,D_2?

% % 2-dim Domain D_1
figure(1)
uu=u(find(SL1+SR1>0),:);
scatter(u(:,1),u(:,2),5,'filled')
hold on
scatter(uu(:,1),uu(:,2),5,'filled')
hold off

% % 2-dim Domain D_2
figure(2)
vv=u(find(SL2+SR2>0),:);
scatter(u(:,1),u(:,2),5,'filled')
hold on
scatter(vv(:,1),vv(:,2),5,'filled')
hold off

% % 3-dim Domain D_1
% uu=u(find(SL1+SR1>0),:);
% scatter3(u(:,1),u(:,2),u(:,3),5,'filled')
% hold on
% scatter3(uu(:,1),uu(:,2),uu(:,3),5,'filled')
% hold off

%% 



% For Figure 7 below, BS model prices and Model 1 prices
%% 2. BS model price (option pricing for SSE 50 ETF) 
% For Figure 7, see "results-50ETF.xls".
clear,clc
MR = xlsread('50ETF-Close.xls'); 
EO = xlsread('ETF-option-9.xls'); 

NumPrc=length(EO);
C_t=zeros(NumPrc,1);
T1=1;

for xx=1:NumPrc

S=MR(xx:243+xx);
X=log(S);    % log price of stock
n=length(S);

% Following is for sigma^2 estimation (quadratic variation mathod)
A=0; 
for i=1:n-1
   A=A+(X(i+1)-X(i))^2;
end
sigma2_hat=A/T1; % sigma^2 estimator


% option pricing below
mu=0; sigma=1; % for standard normal distribution

% Striking prices K
K=2.70; 
% K=2.75; 
% K=2.80; 

r=0.039;
X_0=S(n);
num_days=EO(xx,2);
T=num_days/360; 

sgm=sqrt(sigma2_hat);
d_1=log(S(n)/K)/(sgm*sqrt(T))+(r*sqrt(T))/sgm+0.5*sgm*sqrt(T);
d_2=d_1-sgm*sqrt(T);
Phi_1=cdf('Normal',d_1,mu,sigma);
Phi_2=cdf('Normal',d_2,mu,sigma);
C_t(xx)=S(n)*Phi_1-K*exp(-r*T)*Phi_2;
end
C_t


%% 3. SSE 50 ETF-Estimation and Option Pricing by Model 1(Model 1 Prices)
% For Figure 7, see "results-50ETF.xls".
clear,clc
MR = xlsread('50ETF-Close.xls'); 
EO = xlsread('ETF-option-9.xls'); 

NumPrc=length(EO);
C_1=zeros(length(EO),1);
C_2=zeros(length(EO),1);
C_3=zeros(length(EO),1);

for xx=1:NumPrc

S=MR(xx:243+xx);
X=log(S);    % log price of stock
n=length(S);

R=zeros(n-1,1);

for k=1:n-1
    R(k)=X(k+1)-X(k);   % log return of stock
end

% Step 1. initial parameters, estimate the model by GARCH(1,1) with offset but without risk premia
% First, create a GARCH(1,1) model. 
Mdl_data = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
% Second, fit the GARCH(1,1) model to the data.
EstMdl = estimate(Mdl_data,R,'Display','off');
EstMdl.Offset=0;

NumIterat=5; Index=0; NumPth=50;
alpha_00=zeros(NumPth,1);
alpha_11=zeros(NumPth,1);
beta_11=zeros(NumPth,1);
mu_00=zeros(NumPth,1);
for cnt=1:NumIterat
    
Index=Index+1;
% Step 2. 
[ConVar,Invt] = simulate(EstMdl,numel(R),'NumPaths',NumPth);

for i=1:NumPth
y=R+0.5*ConVar(:,i);
Mdl_y = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
EstMdl_y = estimate(Mdl_y,y,'Display','off');

alpha_00(i)=EstMdl_y.Constant;
alpha_11(i)=EstMdl_y.ARCH{1};
beta_11(i)=EstMdl_y.GARCH{1};
mu_00(i)=EstMdl_y.Offset;
end

aa_0=mean(alpha_00);
aa_1=mean(alpha_11);
bb_1=mean(beta_11);
mu_0=mean(mu_00);
Coff_GARCH=[mu_0 aa_0 aa_1 bb_1];

EstMdl=garch('Constant',aa_0,'GARCH',bb_1,'ARCH',aa_1);
end



% option pricing below
mu=mu_0; alpha_0=aa_0; alpha_1=aa_1; beta=bb_1;
unit_t=1/365; r_1=0.039; r=exp(r_1*unit_t); rho=-(mu-r_1*unit_t);

[ConVar1,Invt1]=Infer_GARCH_11(alpha_0,alpha_1,beta,R,mu);
z=Invt1./sqrt(ConVar1);


t=243+xx; T=t+EO(xx,2);
MC=100000; % numbers of Monte Carlo
% MC=300000; % numbers of Monte Carlo
u=randn(MC,T-t);

x1=zeros(MC,T-t);
x2=zeros(MC,T-t);
ConVF1=zeros(MC,T-t);
ConVF2=zeros(MC,T-t);
for j=1:T-t
    for i=1:MC
    ConVF1(i,j)=ConVolFun(x1(i,1:j),j,alpha_0,alpha_1,beta,ConVar1(n-1),z(n-1));
    x1(i,j)=u(i,j)+rho/ConVF1(i,j)+ConVF1(i,j);
    
    ConVF2(i,j)=ConVolFun(x2(i,1:j),j,alpha_0,alpha_1,beta,ConVar1(n-1),z(n-1));
    x2(i,j)=u(i,j)+rho/ConVF2(i,j);
    end
end

K=2.70;
SL1=sum(ConVF1.*u,2);
SL2=sum(ConVF2.*u,2);
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_1(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

K=2.75;
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_2(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

K=2.80;
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_3(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

end
C_1
C_2
C_3

%% 4. Plot(for Figure 7, see "results-50ETF.xls")
% For Figure 7, see "results-50ETF.xls".
plot(1:length(XX),XX(:,1),'-*',1:length(XX),XX(:,2),'-.+',...
    1:length(XX),XX(:,3),'--k^',...
    'LineWidth',1.5)
legend('Market Prices','BS Prices','Model 1 Prices')
title('K=2.70')
xlim([1 length(XX)])

%%







% Model 2 below, estimation and option pricing for SSE 50 ETF
%% 5. SSE 50 ETF-Estimation and Option Pricing by Model 2
clear,clc
MR = xlsread('50ETF-Close.xls'); 
EO = xlsread('ETF-option-9.xls'); 

C_1=zeros(length(EO),1);
C_2=zeros(length(EO),1);
C_3=zeros(length(EO),1);

for xx=1:length(EO)

S=MR(xx:242+xx); % Rolling Windows
X=log(S);    % log price of stock
n=length(S);

R=zeros(n-1,1);

for k=1:n-1
    R(k)=X(k+1)-X(k);   % log return of stock
end

% Step 1. initial parameters, estimate the model by GARCH(1,1) with offset but without risk premia
% First, create a GARCH(1,1) model. 
Mdl_data = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
% Second, fit the GARCH(1,1) model to the data.
EstMdl = estimate(Mdl_data,R,'Display','off');
EstMdl.Offset=0;

NumIterat=5; Index=0; NumPth=10;
alpha_00=zeros(NumPth,1);
alpha_11=zeros(NumPth,1);
beta_11=zeros(NumPth,1);
mu_00=zeros(NumPth,1);
c_00=zeros(NumPth,1);
for cnt=1:NumIterat
    
Index=Index+1;
% Step 2. 
[ConVar,Invt] = simulate(EstMdl,numel(R),'NumPaths',NumPth);

for i=1:NumPth
    
p = polyfit(ConVar(:,i),R,1);
mu_00(i)=p(2); 
c_00(i)=p(1);
y =R-p(1)*ConVar(:,i)-p(2);
Mdl_y = garch('GARCHLags',1,'ARCHLags',1);
EstMdl_y = estimate(Mdl_y,y,'Display','off');

alpha_00(i)=EstMdl_y.Constant;
alpha_11(i)=EstMdl_y.ARCH{1};
beta_11(i)=EstMdl_y.GARCH{1};

end

aa_0=mean(alpha_00);
aa_1=mean(alpha_11);
bb_1=mean(beta_11);
mu_0=mean(mu_00);
cc_0=mean(c_00);
Coff_GARCH=[mu_0 cc_0 aa_0 aa_1 bb_1];

EstMdl=garch('Constant',aa_0,'GARCH',bb_1,'ARCH',aa_1);
end



% option pricing below
mu=mu_0; c=cc_0; alpha_0=aa_0; alpha_1=aa_1; beta=bb_1;
unit_t=1/365; r_1=0.039; r=exp(r_1*unit_t); 

[ConVar1,Invt1]=Infer_c_GARCH_11(alpha_0,alpha_1,beta,R,mu,c);
z=Invt1./sqrt(ConVar1);


t=242+xx; T=t+EO(xx,2);
MC=100000;
% MC=300000;
u=randn(MC,T-t);

x1=zeros(MC,T-t);
x2=zeros(MC,T-t);
ConVF1=zeros(MC,T-t);
ConVF2=zeros(MC,T-t);
rho1=zeros(MC,T-t);
rho2=zeros(MC,T-t);
for j=1:T-t
    for i=1:MC
    ConVF1(i,j)=ConVolFun(x1(i,1:j),j,alpha_0,alpha_1,beta,ConVar1(n-1),z(n-1));
    rho1(i,j)=-(mu+c*ConVF1(i,j)^2+0.5*ConVF1(i,j)^2-r_1*unit_t);
    x1(i,j)=u(i,j)+rho1(i,j)/ConVF1(i,j)+ConVF1(i,j);
    
    ConVF2(i,j)=ConVolFun(x2(i,1:j),j,alpha_0,alpha_1,beta,ConVar1(n-1),z(n-1));
    rho2(i,j)=-(mu+c*ConVF2(i,j)^2+0.5*ConVF2(i,j)^2-r_1*unit_t);
    x2(i,j)=u(i,j)+rho2(i,j)/ConVF2(i,j);
    end
end

K=2.70;
SL1=sum(ConVF1.*u,2);
SL2=sum(ConVF2.*u,2);
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_1(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

K=2.75;
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_2(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

K=2.80;
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_3(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

end
C_1
C_2
C_3


%% 

