function [ConVar,Invt]=Infer_AR_GARCH_11(alpha_0,alpha_1,beta,R,b)

n=length(R);
ConVar=zeros(n+1,1);
Invt=zeros(n,1);

k=length(b)-2;
g_0=b(1); c=b(end); g=b(2:1:k+1);

ConVar(k+1)=alpha_0;

for i=k+1:n
    Invt(i)=R(i)-g_0-c*ConVar(i)-R(i-k:1:i-1)'*g;
    ConVar(i+1)=alpha_0+alpha_1*Invt(i)*Invt(i)+beta*ConVar(i);
end
ConVar=[ConVar(k+1:2*k,1); ConVar(k+1:n,1)];
Invt=[Invt(k+1:2*k,1); Invt(k+1:n,1)];
end