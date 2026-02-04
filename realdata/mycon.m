function [c,ceq] = mycon(x,p,q)
c = 0;
for j = 1:q
    c = c + pi*x(5*p+j)/(x(5*p+q+j)*x(5*p+2*q+j))^0.5/x(5*p+3*q+j);
end
c = c-1;
ceq = [];