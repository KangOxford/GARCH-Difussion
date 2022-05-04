function ConVF=ConVolFun(x,j,alpha_0,alpha_1,beta,sigma2,z)

% GARCH(1,1) Conditional Volatility Function
if  j<1 
    error('2nd input should be greater than 0');
end

if j==1
   ConVF=sqrt(alpha_0+alpha_1*sigma2*z^2+beta*sigma2);
else
   Q1=0; Q11=1;
   for i=1:j-1
      Q11=Q11*(alpha_1*x(j-i)^2+beta);
      Q1=Q1+Q11;
   end
   Q2=Q11*(alpha_1*sigma2*z^2+beta*sigma2);
   ConVF=sqrt(alpha_0+alpha_0*Q1+Q2);
end

end
