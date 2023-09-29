function [f,g]=rosenbrock(x,a)

n=length(x);
f=0;
g=zeros(n,1);
for k=1:n-1
    T=x(k+1)-x(k)^2;
    S=1-x(k);
    f=f+a*T^2+S^2;
    g(k)=-4*a*x(k)*T-2*S;
    if k>1
        g(k)=g(k)+2*a*(x(k)-x(k-1)^2);
    end
end

g(n)=2*a*(x(n)-x(n-1)^2);

return