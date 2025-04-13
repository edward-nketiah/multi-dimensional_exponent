function [space_time_data,N1,N_b] = fgenerate_data_CRIME(S,T)
%------------------------------------------------------------------------------%
% Modified Fast Algorithm: Based on zhuang, ogata and vere-jones, 2004.        %
% Generate univariate multi-dimensional space-time self-exciting point process.%
%------------------------------------------------------------------------------%
% Generate G0 (G(1), background) based on the thinning methods
rng(1989);
[lamstar,~] = lambda([0,0,0]);
N = poissrnd(lamstar*T*S(1)*S(2));
% x = [rand(N,1)*S(1),rand(N,1)*S(2),rand(N,1)*T];
x = [rand(N,1)*S(1)-S(1)/2,rand(N,1)*S(2)-S(2)/2,rand(N,1)*T];
% x = sortrows(x,3);%There is no need to rearrange here, it does not affect the result: rearrange together at the end
[fun1,~] = lambda(x);
G = x(rand(N,1)<fun1/lamstar,:);
N_b = size(G,1);
G(:,4) = 1;
space_time_data = G;
% [~,fun2x,fun2y,fun2t] = lambda(G);
n = N_b;
% Generate offspring
theta = 0.2;
omega = 0.1;
sigmax = 0.01;
sigmay = 0.1;
while n>0
%---------------------------------------------------------------------------------------------------------%
% Expand the boundary to (-inf, +inf)*(-inf, +inf)*(0, +inf), so that the integral step in the loop
%(corresponding to all the percentage signs in this program) can be omitted%
% Finally, just remove the points beyond the boundary: because its cumulative strength is limited, 
%although more points (outside the boundary) will be calculated, the integral part takes more time. %
% This part needs to be proved.                                                                                          %
lamoff = theta*ones(n,1);
%---------------------------------------------------------------------------------------------------------%
    n = poissrnd(lamoff);
    n1 = sum(n);
    % %Distance from the center point: Actually, fun2x, fun2y and fun2t should be used. 
    %Here is a simplified calculation
    G2 = [randn(n1,2).*[sigmax,sigmay],exprnd(1/omega,n1,1),zeros(n1,1)];
    r = 1;
    G_off = zeros(n1,4);
        while n1>0
            n1 = sum(n>0);
            G_off(r:r+n1-1,1:3) = G2(r:r+n1-1,1:3)+G(n>0,1:3);
            
            n = n-1;               
            r = n1+r;        
        end
%Collection Points
    space_time_data = [space_time_data;G_off];
    G = G_off;
    n = size(G,1);
end
%Remove points outside the boundary
% space_time_data(space_time_data(:,1)<0|space_time_data(:,1)>S(1),:) = [];   
% space_time_data(space_time_data(:,2)<0|space_time_data(:,2)>S(2),:) = []; 
space_time_data(space_time_data(:,1)<-S(1)/2|space_time_data(:,1)>S(1)/2,:) = [];   
space_time_data(space_time_data(:,2)<-S(2)/2|space_time_data(:,2)>S(2)/2,:) = []; 
space_time_data(space_time_data(:,3)>T,:) = []; 
[space_time_data, N1] = sortrows(space_time_data,3);
end

function [fun1,fun2x,fun2y,fun2t] = lambda(x)
%This is a generalized parametrization intensity function 
mu = 5.71;
theta = 0.2;
omega = 0.1;
sigmax = 0.01;
sigmay = 0.1;
fun1 = mu/(2*pi*4.5^2).*exp(-x(:,1).^2/(2*4.5^2)-x(:,2).^2/(2*4.5^2)); 
fun2x = @(s)1/(sqrt(2*pi)*sigmax).*exp(-(s-x(:,1)).^2/(2*sigmax^2));
fun2y = @(s)1/(sqrt(2*pi)*sigmay).*exp(-(s-x(:,2)).^2/(2*sigmay^2));
fun2t = @(s)theta*omega*exp(-omega*(s-x(:,3)));
end