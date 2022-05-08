%% Option Pricing for SPX 
% Underlying Index: SPX

%% Part 1. Estimation for Models


% For section 3.2
%% 1.1. Estimation-Model 1
% For section 3.2
clear,clc
MR = xlsread('SPX-Close.xls'); 
S=MR(4278:4528); % Year 2017/01/03-2017/12/31
% S=MR(4026:4528); % Year 2016/01/04-2017/12/31
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

NumIterat=10; Index=0; NumPth=5;
alpha_00=zeros(NumPth,1);
alpha_11=zeros(NumPth,1);
beta_11=zeros(NumPth,1);
mu_00=zeros(NumPth,1);
for cnt=1:NumIterat
    
Index=Index+1
% Step 2. 
[ConVar,Invt] = simulate(EstMdl,numel(R),'NumPaths',NumPth);

for i=1:NumPth
y=R+0.5*ConVar(:,i);
Mdl_y = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
% EstMdl_y = estimate(Mdl_y,y,'Display','off');
EstMdl_y = estimate(Mdl_y,y,'Display','params');

alpha_00(i)=EstMdl_y.Constant;
alpha_11(i)=EstMdl_y.ARCH{1};
beta_11(i)=EstMdl_y.GARCH{1};
mu_00(i)=EstMdl_y.Offset;
end

aa_0=mean(alpha_00);
aa_1=mean(alpha_11);
bb_1=mean(beta_11);
mu_0=mean(mu_00);
Coff_GARCH=[mu_0 aa_0 aa_1 bb_1]

EstMdl=garch('Constant',aa_0,'GARCH',bb_1,'ARCH',aa_1);
end

%% 1.2 Estimation-Model 2
% For section 3.2
clear,clc
MR = xlsread('SPX-Close.xls'); 
S=MR(4278:4528); % Year 2017/01/03-2017/12/31
% S=MR(4026:4528); % Year 2016/01/04-2017/12/31
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

NumIterat=10; Index=0; NumPth=100;
alpha_00=zeros(NumPth,1);
alpha_11=zeros(NumPth,1);
beta_11=zeros(NumPth,1);
mu_00=zeros(NumPth,1);
c_00=zeros(NumPth,1);
for cnt=1:NumIterat
    
Index=Index+1
% Step 2. 
[ConVar,Invt] = simulate(EstMdl,numel(R),'NumPaths',NumPth);

for i=1:NumPth

p = polyfit(ConVar(:,i),R,1);
mu_00(i)=p(2); 
c_00(i)=p(1);
y =R-p(1)*ConVar(:,i)-p(2);
Mdl_y = garch('GARCHLags',1,'ARCHLags',1);
% EstMdl_y = estimate(Mdl_y,y,'Display','off');
EstMdl_y = estimate(Mdl_y,y,'Display','params');

alpha_00(i)=EstMdl_y.Constant;
alpha_11(i)=EstMdl_y.ARCH{1};
beta_11(i)=EstMdl_y.GARCH{1};

end

aa_0=mean(alpha_00);
aa_1=mean(alpha_11);
bb_1=mean(beta_11);
mu_0=mean(mu_00);
cc_0=mean(c_00);
Coff_GARCH=[mu_0 cc_0 aa_0 aa_1 bb_1]

% StdErr_mu=std(mu_00)/sqrt(length(mu_00))
% StdErr_c=std(c_00)/sqrt(length(c_00))

EstMdl=garch('Constant',aa_0,'GARCH',bb_1,'ARCH',aa_1);
end



%% 1.3 Estimation-Model 3
% For section 3.2
clear,clc
MR = xlsread('SPX-Close.xls'); 
S=MR(4278:4528); % Year 2017/01/03-2017/12/31
% S=MR(4026:4528); % Year 2016/01/04-2017/12/31
X=log(S);    % log price of stock
n=length(S);

R=zeros(n-1,1);
for k=1:n-1
    R(k)=X(k+1)-X(k);   % log return of stock
end

% Step 1. initial parameters, estimate the model by AR(k)-GARCH(1,1) with offset but without risk premia
% First, create a AR(k)-GARCH(1,1) model.
k=1;
Mdl_data = arima('ARLags',k,'Variance',garch(1,1));
% Second, fit the AR(k)-GARCH(1,1) model to the data.
EstMdl = estimate(Mdl_data,R,'Display','off');
EstMdl=EstMdl.Variance;

NumIterat=10; Index=0; NumPth=100;
alpha_00=zeros(NumPth,1);
alpha_11=zeros(NumPth,1);
beta_11=zeros(NumPth,1);
B=zeros(k+2,NumPth);
for cnt=1:NumIterat
    
Index=Index+1
% Step 2. 
[ConVar,Invt] = simulate(EstMdl,numel(R),'NumPaths',NumPth);


for i=1:NumPth
    
X=[zeros(length(R)-k,k) ConVar(k+1:1:end,i)];
for ii=1:k
    X(:,ii)=R(k+1-ii:1:end-ii);
end
Y=R(k+1:1:end);
mdl_lm = fitlm(X,Y);
B(:,i)=mdl_lm.Coefficients.Estimate;
y=Y-[ones(length(R)-k,1) X]*B(:,i);
Mdl_y = garch('GARCHLags',1,'ARCHLags',1);
% EstMdl_y = estimate(Mdl_y,y,'Display','off');
EstMdl_y = estimate(Mdl_y,y,'Display','params');

alpha_00(i)=EstMdl_y.Constant;
alpha_11(i)=EstMdl_y.ARCH{1};
beta_11(i)=EstMdl_y.GARCH{1};

end

aa_0=mean(alpha_00);
aa_1=mean(alpha_11);
bb_1=mean(beta_11);
BB=mean(B,2)'
Coff_GARCH=[aa_0 aa_1 bb_1]

EstMdl=garch('Constant',aa_0,'GARCH',bb_1,'ARCH',aa_1);
end
% StdErr_c=std(B(:,k+2))/sqrt(length(B(:,k+2)))

%% Part 2 Option Pricing 
% Using the parameters estimated from SPX historical data during 2017.
% For section 4, Figure 3.



%% 2.1 Pricing-BS Model-theoretical
% For Figure 3.
clear,clc
MR = xlsread('SPX-Close.xls'); 
S=MR(4278:4528); % Year 2017/01/03-2017/12/31
% S=MR(4026:4528); % Year 2016/01/04-2017/12/31
X=log(S);    % log price of stock
n=length(S);
T1=n/252;

% Following is for sigma^2 estimation (quadratic variation mathod)
A=0;
for i=1:n-1
   A=A+(X(i+1)-X(i))^2;
end
sigma2_hat=A/T1;

% option pricing below
NumPrc=60;
C_t=zeros(NumPrc,1);
K=S(n);  % at-the-money
% K=S(n)*1.05;  % out-of-money
r=0.025;

for xx=1:NumPrc

mu=0; sigma=1; % for standard normal distribution
T=xx/365; 

sgm=sqrt(sigma2_hat);
d_1=log(S(n)/K)/(sgm*sqrt(T))+(r*sqrt(T))/sgm+0.5*sgm*sqrt(T);
d_2=d_1-sgm*sqrt(T);
Phi_1=cdf('Normal',d_1,mu,sigma);
Phi_2=cdf('Normal',d_2,mu,sigma);
C_t(xx)=S(n)*Phi_1-K*exp(-r*T)*Phi_2;

end
C_t


%% 2.2 Pricing-Model 1-theoretical
% For Figure 3.
clear,clc
MR = xlsread('SPX-Close.xls'); 
S=MR(4278:4528); % Year 2017/01/03-2017/12/31
% S=MR(4026:4528); % Year 2016/01/04-2017/12/31
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

NumIterat=10; Index=0; NumPth=5;
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
Coff_GARCH=[mu_0 aa_0 aa_1 bb_1]

EstMdl=garch('Constant',aa_0,'GARCH',bb_1,'ARCH',aa_1);
end
mu=mu_0; alpha_0=aa_0; alpha_1=aa_1; beta=bb_1;

% Option pricing below
NumPrc=60;
C_t=zeros(NumPrc,1);

K=S(n);  % at-the-money
% K=S(n)*1.05;  % out-of-money
unit_t=1/365; r_1=0.025; r=exp(r_1*unit_t); rho=-(mu-r_1*unit_t);
[ConVar1,Invt1]=Infer_GARCH_11(alpha_0,alpha_1,beta,R,mu);
z=Invt1./sqrt(ConVar1);

for xx=1:NumPrc

t=n; T=n+xx;
MC=100000; % Numbers of Monte Carlo
% MC=300000; % Numbers of Monte Carlo
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

SL1=sum(ConVF1.*u,2);
SL2=sum(ConVF2.*u,2);
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_t(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

end
C_t


%% 2.3 Pricing-Model 2-theoretical
% For Figure 3.
clear,clc
MR = xlsread('SPX-Close.xls'); 
S=MR(4278:4528); % Year 2017/01/03-2017/12/31
% S=MR(4026:4528); % Year 2016/01/04-2017/12/31
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

NumIterat=10; Index=0; NumPth=100;
alpha_00=zeros(NumPth,1);
alpha_11=zeros(NumPth,1);
beta_11=zeros(NumPth,1);
mu_00=zeros(NumPth,1);
c_00=zeros(NumPth,1);
for cnt=1:NumIterat
    
Index=Index+1
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
Coff_GARCH=[mu_0 cc_0 aa_0 aa_1 bb_1]

EstMdl=garch('Constant',aa_0,'GARCH',bb_1,'ARCH',aa_1);
end

mu=mu_0; c=cc_0; alpha_0=aa_0; alpha_1=aa_1; beta=bb_1;

% Option pricing below - Model 2
NumPrc=60;
C_t=zeros(NumPrc,1);

K=S(n);  % at-the-money
% K=S(n)*1.25;  % out-of-money
unit_t=1/365; r_1=0.025; r=exp(r_1*unit_t); 
[ConVar1,Invt1]=Infer_c_GARCH_11(alpha_0,alpha_1,beta,R,mu,c);
z=Invt1./sqrt(ConVar1);
for xx=1:NumPrc

t=n; T=n+xx;
MC=100000; % Numbers of Monte Carlo
% MC=300000; % Numbers of Monte Carlo
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

SL1=sum(ConVF1.*u,2);
SL2=sum(ConVF2.*u,2);
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_t(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

end
C_t

% dlmwrite('pricing_Model2.txt',C_t,'delimiter', '\t')

%% 2.4 Pricing-Model 3-theoretical
% For Figure 3.
clear,clc
MR = xlsread('SPX-Close.xls'); 
S=MR(4278:4528); % Year 2017/01/03-2017/12/31
% S=MR(4026:4528); % Year 2016/01/04-2017/12/31
X=log(S);    % log price of stock
n=length(S);

R=zeros(n-1,1);
for k=1:n-1
    R(k)=X(k+1)-X(k);   % log return of stock
end

% Step 1. initial parameters, estimate the model by AR(k)-GARCH(1,1) with offset but without risk premia
% First, create a AR(k)-GARCH(1,1) model.
k=3;
Mdl_data = arima('ARLags',k,'Variance',garch(1,1));
% Second, fit the AR(k)-GARCH(1,1) model to the data.
EstMdl = estimate(Mdl_data,R,'Display','off');
EstMdl=EstMdl.Variance;

NumIterat=10; Index=0; NumPth=100;
alpha_00=zeros(NumPth,1);
alpha_11=zeros(NumPth,1);
beta_11=zeros(NumPth,1);
B=zeros(k+2,NumPth);
for cnt=1:NumIterat
    
Index=Index+1
% Step 2. 
[ConVar,Invt] = simulate(EstMdl,numel(R),'NumPaths',NumPth);

for i=1:NumPth
    
X=[zeros(length(R)-k,k) ConVar(k+1:1:end,i)];
for ii=1:k
    X(:,ii)=R(k+1-ii:1:end-ii);
end
Y=R(k+1:1:end);
mdl_lm = fitlm(X,Y);
B(:,i)=mdl_lm.Coefficients.Estimate;
y=Y-[ones(length(R)-k,1) X]*B(:,i);
Mdl_y = garch('GARCHLags',1,'ARCHLags',1);
EstMdl_y = estimate(Mdl_y,y,'Display','off');

alpha_00(i)=EstMdl_y.Constant;
alpha_11(i)=EstMdl_y.ARCH{1};
beta_11(i)=EstMdl_y.GARCH{1};

end

aa_0=mean(alpha_00);
aa_1=mean(alpha_11);
bb_1=mean(beta_11);
BB=mean(B,2)'
Coff_GARCH=[aa_0 aa_1 bb_1]

EstMdl=garch('Constant',aa_0,'GARCH',bb_1,'ARCH',aa_1);
end

g_0=BB(1); g=BB(2:k+1)'; c=BB(k+2); 
alpha_0=aa_0; alpha_1=aa_1; beta=bb_1;

% Option pricing below - Model 3
NumPrc=60;
C_t=zeros(NumPrc,1);

K=S(n);  % at-the-money
% K=S(n)*1.05;  % out-of-money
unit_t=1/365; r_1=0.025; r=exp(r_1*unit_t); 
[ConVar1,Invt1]=Infer_AR_GARCH_11(alpha_0,alpha_1,beta,R,BB');
z=Invt1./sqrt(ConVar1);
for xx=1:NumPrc

t=n; T=n+xx;
MC=100000; % Numbers of Monte Carlo
% MC=300000; % Numbers of Monte Carlo
u=randn(MC,T-t);

x1=zeros(MC,T-t);
x2=zeros(MC,T-t);
ConVF1=zeros(MC,T-t);
ConVF2=zeros(MC,T-t);
ConRF1=zeros(MC,T-t);
ConRF2=zeros(MC,T-t);
rho1=zeros(MC,T-t);
rho2=zeros(MC,T-t);
for j=1:T-t
    for i=1:MC
    ConVF1(i,j)=ConVolFun(x1(i,1:j),j,alpha_0,alpha_1,beta,ConVar1(n-1),z(n-1));
    ConRF1(i,j)=ConRtnFun(x1(i,1:j),j,BB',R(end-k+1:1:end),[alpha_0;alpha_1;beta],ConVar1(n-1),z(n-1));
    rho1(i,j)=-(ConRF1(i,j)-r_1*unit_t);
    x1(i,j)=u(i,j)+rho1(i,j)/ConVF1(i,j)+ConVF1(i,j);
    
    ConVF2(i,j)=ConVolFun(x2(i,1:j),j,alpha_0,alpha_1,beta,ConVar1(n-1),z(n-1));
    ConRF2(i,j)=ConRtnFun(x2(i,1:j),j,BB',R(end-k+1:1:end),[alpha_0;alpha_1;beta],ConVar1(n-1),z(n-1));
    rho2(i,j)=-(ConRF2(i,j)-r_1*unit_t);
    x2(i,j)=u(i,j)+rho2(i,j)/ConVF2(i,j);
    end
end

SL1=sum(ConVF1.*u,2);
SL2=sum(ConVF2.*u,2);
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_t(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

end
C_t

% dlmwrite('pricing_Model3.txt',C_t,'delimiter', '\t')

%% 2.5 Plot (for Figure 3. See "results-SPX.xls")
% For Figure 3.
plot(1:length(XX),XX(:,1),'-.',1:length(XX),XX(:,5),'-k',...
    1:length(XX),XX(:,6),'--',1:length(XX),XX(:,7),':',...
    'LineWidth',1.5)
legend('B-S','Model 1','Model 2','Model 3')
xlim([1 length(XX)])


%%


%% 2.6 The effects of return rate mu on option price (An informal analysis)

% Option pricing below
K=S(n);  % at-the-money
% K=S(n)*1.05;  % out-of-money
unit_t=1/365; r_1=0.025; r=exp(r_1*unit_t); 

t=n; T=n+30;
MC=1000;  % Numbers of Monte Carlo
% MC=10000;  % Numbers of Monte Carlo 
u=randn(MC,T-t);

x1=zeros(MC,T-t);
x2=zeros(MC,T-t);
ConVF1=zeros(MC,T-t);
ConVF2=zeros(MC,T-t);
lgmu=-6:0.05:0;  % the return rate mu from 10^(-6) to 10^(0)=1. The empirical reurn rate mu=10^(-4).
mu=10.^(lgmu)';
C_t=zeros(length(mu),1);

for xx=1:length(mu)

rho=-(mu(xx)-r_1*unit_t);
[ConVar1,Invt1]=Infer_GARCH_11(alpha_0,alpha_1,beta,R,mu(xx));
z=Invt1./sqrt(ConVar1);
for j=1:T-t
    for i=1:MC
    ConVF1(i,j)=ConVolFun(x1(i,1:j),j,alpha_0,alpha_1,beta,ConVar1(n-1),z(n-1));
    x1(i,j)=u(i,j)+rho/ConVF1(i,j)+ConVF1(i,j);
    
    ConVF2(i,j)=ConVolFun(x2(i,1:j),j,alpha_0,alpha_1,beta,ConVar1(n-1),z(n-1));
    x2(i,j)=u(i,j)+rho/ConVF2(i,j);
    end
end

SL1=sum(ConVF1.*u,2);
SL2=sum(ConVF2.*u,2);
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_t(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;
[xx C_t(xx)]
end
C_t
plot(lgmu,C_t,'-k','LineWidth',1)



%% Part 3. Pricing by rolling windows, camparing with market prices, empirical applications

% For section 5.1, Figure 4 and Figure 5.




%% 3.1 BS model price-option pricing C2700
% For Figure 5
clear,clc
MR = xlsread('SPX-Close.xls'); 
EO = xlsread('SPX-option-C2700-0316.xls'); 

NumPrc=length(EO);
C_t=zeros(NumPrc,1);
T1=1;

for xx=1:NumPrc

S=MR(4281+xx:4533+xx);
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
K=2700; 
r=0.025;
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


%% 3.2 SPX-Estimation and Option Pricing by Model 1; C2700
% For Figure 5
clear,clc
MR = xlsread('SPX-Close.xls'); 
EO = xlsread('SPX-option-C2700-0316.xls'); 

NumPrc=length(EO);
C_t=zeros(NumPrc,1);

for xx=1:NumPrc

S=MR(4281+xx:4533+xx);
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

NumIterat=10; Index=0; NumPth=20;
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
unit_t=1/365; r_1=0.025; r=exp(r_1*unit_t); rho=-(mu-r_1*unit_t);

[ConVar1,Invt1]=Infer_GARCH_11(alpha_0,alpha_1,beta,R,mu);
z=Invt1./sqrt(ConVar1);


t=242+xx; T=t+EO(xx,2);
MC=100000; % Numbers of Monte Carlo
% MC=200000; % Numbers of Monte Carlo
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

K=2700;
SL1=sum(ConVF1.*u,2);
SL2=sum(ConVF2.*u,2);
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_t(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

end
C_t

%% 3.3 plot, Market-BS-Model 1 (See "results-SPX-option-C2700-0316.xls")
% For Figure 5.
plot(1:length(XX),XX(:,1),'-*',1:length(XX),XX(:,2),'--x',...
    1:length(XX),XX(:,3),'-.k+',...
    'LineWidth',1.5,'MarkerSize',7)
legend('Market Prices','BS Prices','Model 1 Prices')
title('K=2700')
xlim([1 length(XX)])


%%


%% 3.4 BS model price-option pricing C2650.
% For Figure 4.
clear,clc
MR = xlsread('SPX-Close.xls'); 
EO = xlsread('SPX-option-C2650-0216.xls'); 

NumPrc=length(EO);
C_t=zeros(NumPrc,1);
T1=1;

for xx=1:NumPrc

S=MR(4257+xx:4508+xx);
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
K=2650; 
r=0.025;
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

%% 3.5 SPX-Estimation and Option Pricing by Model 1. (C2650)
% For Figure 4.
clear,clc
MR = xlsread('SPX-Close.xls'); 
EO = xlsread('SPX-option-C2650-0216.xls'); 

NumPrc=length(EO);
C_t=zeros(NumPrc,1);

for xx=1:NumPrc

S=MR(4257+xx:4508+xx);
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

NumIterat=10; Index=0; NumPth=20;
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
unit_t=1/365; r_1=0.025; r=exp(r_1*unit_t); rho=-(mu-r_1*unit_t);

[ConVar1,Invt1]=Infer_GARCH_11(alpha_0,alpha_1,beta,R,mu);
z=Invt1./sqrt(ConVar1);


t=242+xx; T=t+EO(xx,2);
MC=100000; % Numbers of Monte Carlo
% MC=300000; % Numbers of Monte Carlo
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

K=2650;
SL1=sum(ConVF1.*u,2);
SL2=sum(ConVF2.*u,2);
SR1=log(S(n)/K)+r_1*unit_t*(T-t)+0.5*sum(ConVF1.^2,2);
SR2=log(S(n)/K)+r_1*unit_t*(T-t)-0.5*sum(ConVF2.^2,2);
count1=length(find(SL1+SR1>0));
count2=length(find(SL2+SR2>0));
Phi_D1=count1/MC;
Phi_D2=count2/MC;
C_t(xx)=S(n)*Phi_D1-K*r^(t-T)*Phi_D2;

end
C_t

%% 3.6 plot, Market-BS-Model 1 (See "results-SPX-option-C2650-0216.xls")
% For Figure 4.
plot(1:length(XX),XX(:,1),'-*',1:length(XX),XX(:,2),'--x',...
    1:length(XX),XX(:,3),'-.k+',...
    'LineWidth',1.5,'MarkerSize',7)
legend('Market Prices','BS Prices','Model 1 Prices')
title('K=2650')
xlim([1 length(XX)])

%%



