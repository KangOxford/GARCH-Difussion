function ConRtn=ConRtnFun(x,j,b,R,alpha,sigma2,z)

if  j<1 
    error('2nd input should be greater than 0');
end

k=length(b)-2;
g_0=b(1); c=b(end); g=b(2:1:k+1);

ConVF=zeros(j,1);
for i=1:j
    ConVF(i)=ConVolFun(x(1:i),i,alpha(1),alpha(2),alpha(3),sigma2,z);
end

mu=zeros(j,1);
RR=zeros(j,1);
mu(1)=g_0+(c+0.5)*ConVF(1)^2+R(end-k+1:1:end)'*g(end:-1:1);

if j>1 && j<=k
    for p=1:j-1
        RR(p)=mu(p)-0.5*ConVF(p)^2+ConVF(p)*x(p);
        mu(p+1)=g_0+(c+0.5)*ConVF(p+1)^2+[R(end-k+p+1:1:end);RR(1:p)]'*g(end:-1:1);
    end
elseif j>=k+1
    for p=1:k-1
        RR(p)=mu(p)-0.5*ConVF(p)^2+ConVF(p)*x(p);
        mu(p+1)=g_0+(c+0.5)*ConVF(p+1)^2+[R(end-k+p+1:1:end);RR(1:p)]'*g(end:-1:1);
    end
    for p=k:j-1
        RR(p)=mu(p)-0.5*ConVF(p)^2+ConVF(p)*x(p);
        mu(p+1)=g_0+(c+0.5)*ConVF(p+1)^2+RR(p-k+1:p)'*g(end:-1:1);
    end
end
    
ConRtn=mu(j);

end

