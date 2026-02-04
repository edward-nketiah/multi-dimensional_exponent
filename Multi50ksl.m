
S = load('x_mles_50k_data.mat'); 

x_mles_kept = S.x_mles_kept;
p = 0; 
q = 1;

% model sizes
p_tr = p+1; 
q_tr = q+1;
n_params = 5*(p+1) + 4*(q+1); 

% Raw mean/var 
x_mean = mean(x_mles_kept, 1);
x_var  = var( x_mles_kept, 0, 1);

fprintf('\n================ estimates =================\n');
for k = 1:numel(x_mean)
    fprintf('param %2d: mean = %+ .6f, var = %.6e\n', k, x_mean(k), x_var(k));
end
fprintf('========================================================\n');

% Transformed arrays sized to N_kept 
n_kept = size(x_mles_kept, 1);

mu_bar1_all=NaN(n_kept,1); 
mu_bar2_all=NaN(n_kept,1);
c_all=NaN(n_kept,1); 
d_all=NaN(n_kept,1);
sigma_mu1_all=NaN(n_kept,1); 
sigma_mu2_all=NaN(n_kept,1);

theta1_all=NaN(n_kept,1); 
omega1_all=NaN(n_kept,1); 
sigma_x1_all=NaN(n_kept,1); 
sigma_y1_all=NaN(n_kept,1);
theta2_all=NaN(n_kept,1); 
omega2_all=NaN(n_kept,1); 
sigma_x2_all=NaN(n_kept,1); 
sigma_y2_all=NaN(n_kept,1);

for ii = 1:n_kept
    xs = x_mles_kept(ii,:);

    % Background
    t1b=(2*xs(2*p_tr+1))^(-0.5); % idx 3
    t2b=(2*xs(4*p_tr+1))^(-0.5); % idx 5
    mu_bar1_all(ii)=2*pi*t1b*t1b*xs(1);
    mu_bar2_all(ii)=2*pi*t2b*t2b*xs(1);
    c_all(ii)=xs(p_tr+1);        % idx 2
    d_all(ii)=xs(3*p_tr+1);      % idx 4
    sigma_mu1_all(ii)=t1b; 
    sigma_mu2_all(ii)=t2b;

    % Triggering (comp 1)
    t1g1=(2*xs(5*p_tr+q_tr+1))^(-0.5);   % idx 8
    t2g1=(2*xs(5*p_tr+2*q_tr+1))^(-0.5); % idx 10
    theta1_all(ii)=xs(5*p_tr+1)*2*pi*t1g1*t2g1/xs(5*p_tr+3*q_tr+1); % 6 / 12
    omega1_all(ii)=xs(5*p_tr+3*q_tr+1);  % 12
    sigma_x1_all(ii)=t1g1; sigma_y1_all(ii)=t2g1;

    % Triggering (comp 2)
    t1g2=(2*xs(5*p_tr+q_tr+2))^(-0.5);   % idx 9
    t2g2=(2*xs(5*p_tr+2*q_tr+2))^(-0.5); % idx 11
    theta2_all(ii)=xs(5*p_tr+2)*2*pi*t1g2*t2g2/xs(5*p_tr+3*q_tr+2); % 7 / 13
    omega2_all(ii)=xs(5*p_tr+3*q_tr+2);  % 13
    sigma_x2_all(ii)=t1g2; 
    sigma_y2_all(ii)=t2g2;
end

% Means/vars
mu_bar1_mean = mean(mu_bar1_all); 
mu_bar1_var = var(mu_bar1_all,0);
mu_bar2_mean = mean(mu_bar2_all); 
mu_bar2_var = var(mu_bar2_all,0);
c_mean = mean(c_all);             
c_var = var(c_all,0);
d_mean = mean(d_all);             
d_var = var(d_all,0);
sigma_mu1_mean = mean(sigma_mu1_all); 
sigma_mu1_var = var(sigma_mu1_all,0);
sigma_mu2_mean = mean(sigma_mu2_all); 
sigma_mu2_var = var(sigma_mu2_all,0);

theta1_mean = mean(theta1_all);   
theta1_var = var(theta1_all,0);
omega1_mean = mean(omega1_all);   
omega1_var = var(omega1_all,0);
sigma_x1_mean = mean(sigma_x1_all); 
sigma_x1_var = var(sigma_x1_all,0);
sigma_y1_mean = mean(sigma_y1_all); 
sigma_y1_var = var(sigma_y1_all,0);

theta2_mean = mean(theta2_all);   
theta2_var = var(theta2_all,0);
omega2_mean = mean(omega2_all);   
omega2_var = var(omega2_all,0);
sigma_x2_mean = mean(sigma_x2_all); 
sigma_x2_var = var(sigma_x2_all,0);
sigma_y2_mean = mean(sigma_y2_all); 
sigma_y2_var = var(sigma_y2_all,0);

fprintf('\n=========== Transformed estimates ===========\n', n_kept);
fprintf('Background:\n');
fprintf('  mu_bar1: mean = %.6f, var = %.6e\n', mu_bar1_mean, mu_bar1_var);
fprintf('  mu_bar2: mean = %.6f, var = %.6e\n', mu_bar2_mean, mu_bar2_var);
fprintf('  c      : mean = %.6f, var = %.6e\n', c_mean,       c_var);
fprintf('  d      : mean = %.6f, var = %.6e\n', d_mean,       d_var);
fprintf('  sigma_mu1: mean = %.6f, var = %.6e\n', sigma_mu1_mean, sigma_mu1_var);
fprintf('  sigma_mu2: mean = %.6f, var = %.6e\n', sigma_mu2_mean, sigma_mu2_var);

fprintf('Triggering (component 1):\n');
fprintf('  theta1 : mean = %.6f, var = %.6e\n', theta1_mean, theta1_var);
fprintf('  omega1 : mean = %.6f, var = %.6e\n', omega1_mean, omega1_var);
fprintf('  sigma_x1: mean = %.6f, var = %.6e\n', sigma_x1_mean, sigma_x1_var);
fprintf('  sigma_y1: mean = %.6f, var = %.6e\n', sigma_y1_mean, sigma_y1_var);

fprintf('Triggering (component 2):\n');
fprintf('  theta2 : mean = %.6f, var = %.6e\n', theta2_mean, theta2_var);
fprintf('  omega2 : mean = %.6f, var = %.6e\n', omega2_mean, omega2_var);
fprintf('  sigma_x2: mean = %.6f, var = %.6e\n', sigma_x2_mean, sigma_x2_var);
fprintf('  sigma_y2: mean = %.6f, var = %.6e\n', sigma_y2_mean, sigma_y2_var);
fprintf('=======================================================================\n');



