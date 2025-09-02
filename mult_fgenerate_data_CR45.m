function [space_time_data, N1, N_b] = mult_fgenerate_data_CR45(S, T)
%------------------------------------------------------------------------------%
% Multi-component recursive self-exciting spatiotemporal point process model   %
% Two-component offspring generation using double while-loop structure         %
%------------------------------------------------------------------------------%

rng(1985); %  Reproducibilityfor 5k-1985 and 1981 for 1350

% Step 1: Generate background events
lamstar = lambda([0,0]);
N = poissrnd(lamstar * T * S(1) * S(2));
x = [rand(N,1)*S(1) - S(1)/2, rand(N,1)*S(2) - S(2)/2, rand(N,1)*T];

fun1 = lambda(x);
G = x(rand(N,1) < fun1 / lamstar, :);
N_b = size(G,1);
G(:,4) = 1;  % background label
space_time_data = G;

n = N_b;

% ---------------- Component 1 parameters ----------------
theta1 = 0.2; 
omega1 = 0.1;
sigmax1 = 0.01; 
sigmay1 = 0.1; 
rho1 = 0.0;
Sigma1 = [sigmax1^2, rho1*sigmax1*sigmay1;
          rho1*sigmax1*sigmay1, sigmay1^2];
L1 = chol(Sigma1, 'lower');

% ---------------- Component 2 parameters ----------------
theta2 = 0.3; 
omega2 = 0.2;
sigmax2 = 0.2; 
sigmay2 = 0.05; 
rho2 = 0.0;
Sigma2 = [sigmax2^2, rho2*sigmax2*sigmay2;
          rho2*sigmax2*sigmay2, sigmay2^2];
L2 = chol(Sigma2, 'lower');

% Recursive generation
while n > 0
    % Total offspring per parent for each component
    lamoff = (theta1+theta2) * ones(n,1);
    n = poissrnd(lamoff);
    n1 = sum(n);

    % Prepare correlated spatial noise
    Z1 = randn(n1, 2);
    G21 = [Z1 * L1', exprnd(1/omega1, n1, 1), zeros(n1, 1)];
    Z2 = randn(n1, 2);
    G22 = [Z2 * L2', exprnd(1/omega2, n1, 1), zeros(n1, 1)];

    % Offspring container
    G_off = zeros(n1, 4);
    U = rand(n1,1);
    temp = U < theta1/(theta1+theta2);
    temp1 = ~temp;
    r = 1;
    while n1 > 0
        n1 = sum(n > 0);
        if n1 == 0
            break;
        end

        % Generate offspring by adding noise to parent events
        G_off(r:r+n1-1, 1:3) = temp(r:r+n1-1).*G21(r:r+n1-1, 1:3)...
            +temp1(r:r+n1-1).*G22(r:r+n1-1, 1:3) + G(n > 0, 1:3);
        r = r + n1;
        n = n - 1;
    end

    % Update dataset and parent pool
    space_time_data = [space_time_data; G_off];
    G = G_off;
    n = size(G, 1);
end

% Step 3: Remove out-of-bound events
space_time_data(space_time_data(:,1) < -S(1)/2 | space_time_data(:,1) > S(1)/2, :) = [];
space_time_data(space_time_data(:,2) < -S(2)/2 | space_time_data(:,2) > S(2)/2, :) = [];
space_time_data(space_time_data(:,3) > T, :) = [];

% Step 4: Sort by time
[space_time_data, N1] = sortrows(space_time_data, 3);
end

% ---------------- Lambda: Background intensity ----------------
function fun1 = lambda(x)
mu = 5.71;
rho = 0.0;
sigma = 4.5;

fun1 = mu / (2*pi*sigma^2) .* exp(1/(1 - rho^2) * ...
    (-x(:,1).^2 / (2*sigma^2) - x(:,2).^2 / (2*sigma^2) + ...
     rho * x(:,1).*x(:,2) / sigma^2));
end
