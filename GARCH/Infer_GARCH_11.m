function [ConVar,Invt]=Infer_GARCH_11(alpha_0,alpha_1,beta,R,mu)

n=length(R);
ConVar=zeros(n+1,1);
Invt=zeros(n,1);

ConVar(1)=alpha_0;
for i=1:n
    Invt(i)=R(i)-mu+0.5*ConVar(i);
    ConVar(i+1)=alpha_0+alpha_1*Invt(i)*Invt(i)+beta*ConVar(i);
end
ConVar=ConVar(1:n,1);

end