function [data,lambda,k] = di_method(a,x1,y1,t1,X,Y,T,S,p,q)
x_1 = S(1);%Rectangular boundary, panned to the origin 
x_2 = S(2);
y_1 = S(3);
y_2 = S(4);
t_1 = t1(1);
t_2 = t1(end);
n = length(x1);
mu0 = @(b1,b2,b3,b4)exp(-b2*(x1-b1).^2-b4*(y1-b3).^2);
mu = zeros(p,n);
for j = 1:p
    mu(j,:) = mu0(a(p+j),a(2*p+j),a(3*p+j),a(4*p+j));
end
g0 = @(b1,b2,b3,x,y,t)exp(-b1*x.^2-b2*y.^2-b3*t);
G = zeros(q,n);
for j = 1:q
    G(j,2:end) = sum(g0(a(5*p+q+j),a(5*p+2*q+j),a(5*p+3*q+j),X,Y,T));
end           
%h0 = @(b)(exp(-b*t_2)-exp(-b*t_1))/b;
A = zeros(p,1);
for j = 1:p
    A(j) = (normcdf(x_2,a(p+j),(2*a(2*p+j)).^-0.5)-normcdf(x_1,a(p+j),(2*a(2*p+j)).^-0.5))...
        *(normcdf(y_2,a(3*p+j),(2*a(4*p+j)).^-0.5)-normcdf(y_1,a(3*p+j),(2*a(4*p+j)).^-0.5))...
        *(t_2-t_1)*pi*(a(2*p+j)*a(4*p+j))^-0.5;
end
h = @(b,t)(1-exp(-b*(t_2-t)))/b;
H = zeros(q,n);
B = zeros(q,n);
for j = 1:q
    H(j,:) = h(a(5*p+3*q+j),t1);
    B(j,:) = (normcdf(x_2,x1,(2*a(5*p+q+j)).^-0.5)-normcdf(x_1,x1,(2*a(5*p+q+j)).^-0.5))...
        .*(normcdf(y_2,y1,(2*a(5*p+2*q+j)).^-0.5)-normcdf(y_1,y1,(2*a(5*p+2*q+j)).^-0.5))...
*pi*(a(5*p+q+j)*a(5*p+2*q+j))^-0.5;
end
% k = (-a(1:p)*A+sum(a(6*p+1:6*p+q)*(B.*H)))/S(2)/S(4);%calculate k
lambda = a(1:p)*mu + a(5*p+1:5*p+q)*G;

k = 329.788/(8*15*364.9);

temp = (200*lambda.^-1./sum(lambda.^-1))>rand(1,n);

data = [x1(temp),y1(temp),t1(temp)];
