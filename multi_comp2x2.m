clearvars
clc
% dbstop if error
S = [20, 20];
T = [0, 5000];  % At 1350, there are about 5,000 samples
[space_time_data, N1, N_b] = mult_fgenerate_data_CR45(S, T(2));

% Trim edges to avoid boundary effects
Bound = size(space_time_data, 1);
N1 = N1(2001:Bound-2000);
N_b1 = sum(N1 <= N_b);
space_time_data = space_time_data(2001:Bound-2000, :);
space_time_data(:, 1:2) = space_time_data(:, 1:2) + 10;
space_time_data = space_time_data(:, 1:3);

% Extract x, y, t and shift time
x1 = space_time_data(:, 1);
y1 = space_time_data(:, 2);
t1 = space_time_data(:, 3);
t1 = t1 - t1(1);  % Time shift

n = length(x1);
n1 = 500;

% Create lagged spatial-temporal matrices
X = zeros(n1-1, n-1);
Y = zeros(n1-1, n-1);
Tmat = zeros(n1-1, n-1);
for j = 1:n1-1
    X(j, :) = x1(2:end) - [-Inf*ones(j-1, 1); x1(1:end-j)];
    Y(j, :) = y1(2:end) - [-Inf*ones(j-1, 1); y1(1:end-j)];
    Tmat(j, :) = t1(2:end) - [-Inf*ones(j-1, 1); t1(1:end-j)];
end

% Optimization setup
p = 0+0;
q = 0+1;
%rng(888);
int_a = [rand(1, 5*(p+1)), rand(1, 4*(q+1))];

lb = zeros(1, 5*(p+1) + 4*(q+1));

[x, fval] = fmincon(@(a) li_fun1(a, x1, y1, t1, X, Y, Tmat, [0 20 0 20], p+1, q+1), ...
                        int_a, [], [], [], [], lb, [], @(x) mycon(x, p+1, q+1));


% --------------------------
% Plotting
% --------------------------

mu0 = @(b0, b1, b2, b3, b4, u, v) b0 .* exp(-b2.*(u-b1).^2 - b4.*(v-b3).^2);
mu_cor = -10:0.001:30;
mu1_cor = -10:0.001:30;
n = length(mu_cor);
n1 = length(mu1_cor);

p = p+1;
q = q+1;

mu_x = zeros(p, n);
mu_y = zeros(p, n1);
nor_mu = 0;

for j = 1:p
    mu_x(j, :) = mu0(x(j), x(p+j), x(2*p+j), x(3*p+j), x(4*p+j), mu_cor, x(3*p+j)) * sqrt(pi/x(4*p+j));
    mu_y(j, :) = mu0(x(j), x(p+j), x(2*p+j), x(3*p+j), x(4*p+j), x(p+j), mu1_cor) * sqrt(pi/x(2*p+j));
    nor_mu = nor_mu + pi * x(j) / sqrt(x(2*p+j) * x(4*p+j));
end

% g function setup
g0 = @(b0, b1, b2, b3, u, v, s) b0 .* exp(-b1*u.^2 - b2*v.^2 - b3*s);
g_cor = 0:0.001:50;
g1_cor = -0.5:0.001:0.5;
g2_cor = -0.5:0.0001:0.5;

n = length(g_cor);
n1 = length(g1_cor);
n2 = length(g2_cor);

g_x = zeros(q, n1);
g_y = zeros(q, n2);
g_t = zeros(q, n);
nor_g = 0;

for j = 1:q
    g_x(j, :) = g0(x(5*p+j), x(5*p+q+j), x(5*p+2*q+j), x(5*p+3*q+j), g1_cor, 0, 0) ...
                * sqrt(pi/x(5*p+2*q+j)) / x(5*p+3*q+j);
    g_y(j, :) = g0(x(5*p+j), x(5*p+q+j), x(5*p+2*q+j), x(5*p+3*q+j), 0, g2_cor, 0) ...
                * sqrt(pi/x(5*p+q+j)) / x(5*p+3*q+j);
    g_t(j, :) = g0(x(5*p+j), x(5*p+q+j), x(5*p+2*q+j), x(5*p+3*q+j), 0, 0, g_cor) ...
                * pi / sqrt(x(5*p+q+j) * x(5*p+2*q+j));
    nor_g = nor_g + pi * x(5*p+j) / sqrt(x(5*p+q+j) * x(5*p+2*q+j)) / x(5*p+3*q+j);
end

% Plot results
figure(1)
subplot(2,3,1); 
plot(mu_cor, sum(mu_x, 1)/nor_mu); 
title('mu_x')
subplot(2,3,2); 
plot(mu1_cor, sum(mu_y, 1)/nor_mu); 
title('mu_y')
subplot(2,3,4); 
plot(g1_cor, sum(g_x, 1)/nor_g); 
title('g_x')
subplot(2,3,5); 
plot(g2_cor, sum(g_y, 1)/nor_g); 
title('g_y')
subplot(2,3,6); 
plot(g_cor, sum(g_t, 1)/nor_g); 
title('g_t')

%Show the results
% temp1 = (2*x(2*p+1))^(-0.5);
% temp2 = (2*x(4*p+1))^(-0.5);
% fprintf('\n\n----True Parameter of background function----\n')
% fprintf('mu_bar = %.4f, mu_x1 = %.4f, mu_y1 = %.4f, sigma_x1 = %.4f, sigma_y1 = %.4f',...
%     5.71, 10, 10, 4.5, 4.5);
% fprintf('\n\n----Parameter estimation of background function----\n')
% fprintf('mu_bar = %.4f, mu_x1 = %.4f, mu_y1 = %.4f, sigma_x1 = %.4f, sigma_y1 = %.4f',...
%     2*pi*temp1*temp2*x(1), x(p+1), x(3*p+1), temp1, temp2);
% 
% temp1 = (2*x(5*p+q+1))^(-0.5);
% temp2 = (2*x(5*p+2*q+1))^(-0.5);
% temp3 = (2*x(5*p+q+2))^(-0.5);
% temp4 = (2*x(5*p+2*q+2))^(-0.5);
% temp5 = x(5*p+1)*2*pi*temp1*temp2/x(5*p+3*q+1);
% temp6 = x(5*p+2)*2*pi*temp3*temp4/x(5*p+3*q+2);
% fprintf('\n\n----True Parameter of triggering function----\n')
% fprintf('theta1 = %.4f, omega1 = %.4f, sigma_x1 = %.4f, sigma_y1 = %.4f\n',...
%     0.2, 0.1, 0.01, 0.1);
% fprintf('theta2 = %.4f, omega2 = %.4f, sigma_x2 = %.4f, sigma_y2 = %.4f\n',...
%     0.3, 0.2, 0.2, 0.05);
% fprintf('\n\n----Parameter estimation of triggering function----\n')
% fprintf('theta1 = %.4f, omega1 = %.4f, sigma_x1 = %.4f, sigma_y1 = %.4f\n',...
%     temp5, x(5*p+3*q+1), temp1, temp2);
% fprintf('theta2 = %.4f, omega2 = %.4f, sigma_x2 = %.4f, sigma_y2 = %.4f\n',...
%     temp6, x(5*p+3*q+2), temp3, temp4);
% fprintf('-------------------------------------------------------------------\n')


temp1 = (2*x(2*p+1))^(-0.5);
temp2 = (2*x(4*p+1))^(-0.5);
fprintf('\n\n----True Parameter of background function----\n')
fprintf('mu_bar = %.4f, c = %.4f, d = %.4f, sigma_mu1 = %.4f, sigma_mu2 = %.4f',...
    5.71, 10, 10, 4.5, 4.5);
fprintf('\n\n----Parameter estimation of background function----\n')
fprintf('mu_bar1 = %.4f,mu_bar2 = %.4f, c = %.4f, d = %.4f, sigma_mu1 = %.4f, sigma_mu2 = %.4f',...
    2*pi*temp1*temp1*x(1), 2*pi*temp2*temp2*x(1), x(p+1), x(3*p+1), temp1, temp2);

temp1 = (2*x(5*p+q+1))^(-0.5);
temp2 = (2*x(5*p+2*q+1))^(-0.5);
temp3 = (2*x(5*p+q+2))^(-0.5);
temp4 = (2*x(5*p+2*q+2))^(-0.5);
temp5 = x(5*p+1)*2*pi*temp1*temp2/x(5*p+3*q+1);
temp6 = x(5*p+2)*2*pi*temp3*temp4/x(5*p+3*q+2);
fprintf('\n\n----True Parameter of triggering function----\n')
fprintf('theta1 = %.4f, omega1 = %.4f, sigma_x1 = %.4f, sigma_y1 = %.4f\n',...
    0.2, 0.1, 0.01, 0.1);
fprintf('theta2 = %.4f, omega2 = %.4f, sigma_x2 = %.4f, sigma_y2 = %.4f\n',...
    0.3, 0.2, 0.2, 0.05);
fprintf('\n\n----Parameter estimation of triggering function----\n')
fprintf('theta1 = %.4f, omega1 = %.4f, sigma_x1 = %.4f, sigma_y1 = %.4f\n',...
    temp5, x(5*p+3*q+1), temp1, temp2);
fprintf('theta2 = %.4f, omega2 = %.4f, sigma_x2 = %.4f, sigma_y2 = %.4f\n',...
    temp6, x(5*p+3*q+2), temp3, temp4);
fprintf('-------------------------------------------------------------------\n')


