function [space_time_data, N1, N_b] = fgenerate_data_CRIME21(S, T)
%------------------------------------------------------------------------------%
% Modified Fast Algorithm: Based on Zhuang, Ogata and Vere-Jones (2004).       %
% Generate univariate multi-dimensional space-time self-exciting point process.%
%------------------------------------------------------------------------------%

% Generate G0 (background) using thinning
rng(1989); %0-1982; 1987 for 1350 and 1989 for 15500
lamstar = lambda([0,0]);
N = poissrnd(lamstar * T * S(1) * S(2));
x = [rand(N,1)*S(1) - S(1)/2, rand(N,1)*S(2) - S(2)/2, rand(N,1)*T];

fun1 = lambda(x);
G = x(rand(N,1) < fun1 / lamstar, :);
N_b = size(G, 1);
G(:,4) = 1;
space_time_data = G;

n = N_b;

% Parameters
theta = 0.2;
omega = 0.1;
sigmax = 0.01;
sigmay = 0.1;
rho = 0.4; %0.01

% Covariance matrix
Sigma = [sigmax^2, rho * sigmax * sigmay;
         rho * sigmax * sigmay, sigmay^2];
L = chol(Sigma,'lower');
%
while n > 0
    lamoff = theta * ones(n,1);
    n = poissrnd(lamoff);
    n1 = sum(n);

    % Prepare correlated spatial noise
    Z = randn(n1, 2);        % independent normal samples
    G2 = [Z * L', exprnd(1/omega, n1, 1), zeros(n1, 1)];

    % Offspring container
    G_off = zeros(n1, 4);
    r = 1;

    while n1 > 0
        n1 = sum(n > 0);
        if n1 == 0
            break;
        end

        % Generate offspring by adding noise to parent events
        G_off(r:r+n1-1, 1:3) = G2(r:r+n1-1, 1:3) + G(n > 0, 1:3);
        r = r + n1;
        n = n - 1;
    end

    % Append offspring to dataset
    space_time_data = [space_time_data; G_off];
    G = G_off;
    n = size(G, 1);
end

% Remove points outside boundaries
space_time_data(space_time_data(:,1) < -S(1)/2 | space_time_data(:,1) > S(1)/2, :) = [];
space_time_data(space_time_data(:,2) < -S(2)/2 | space_time_data(:,2) > S(2)/2, :) = [];
space_time_data(space_time_data(:,3) > T, :) = [];

% Sort by time
[space_time_data, N1] = sortrows(space_time_data, 3);
end

function fun1 = lambda(x)
%This is a generalized parametrization intensity function 
mu = 5.71;
rho = 0.4;

fun1 = mu/(2*pi*4.5^2).*exp(1/(1-rho^2)*(-x(:,1).^2/(2*4.5^2)-x(:,2).^2/(2*4.5^2)...
    +rho*x(:,1).*x(:,2)/(4.5*4.5))); 

end