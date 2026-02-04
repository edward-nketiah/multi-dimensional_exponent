S = load('x_mles_5k_singledata.mat');   

x_mles_kept = S.x_mles_kept;
p = 0; 
q = 0;

% model sizes
p_tr = p+1; 
q_tr = q+1;
%n_params = 5*(p+1) + 4*(q+1); 

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

mu_bar1_all   = NaN(n_kept,1);
mu_bar2_all   = NaN(n_kept,1);
c_all         = NaN(n_kept,1);
d_all         = NaN(n_kept,1);
sigma_mu1_all = NaN(n_kept,1);
sigma_mu2_all = NaN(n_kept,1);

theta_all     = NaN(n_kept,1);
omega_all     = NaN(n_kept,1);
sigma_x_all   = NaN(n_kept,1);
sigma_y_all   = NaN(n_kept,1);

for ii = 1:n_kept
    xs = x_mles_kept(ii,:);
    if any(~isfinite(xs)) 
        continue 
    end  
    % Background 
    t1b = (2*xs(2*p_tr+1))^(-0.5);
    t2b = (2*xs(4*p_tr+1))^(-0.5);

    mu_bar1_all(ii)   = 2*pi*t1b*t1b*xs(1);
    mu_bar2_all(ii)   = 2*pi*t2b*t2b*xs(1);
    c_all(ii)         = xs(p_tr+1);
    d_all(ii)         = xs(3*p_tr+1);
    sigma_mu1_all(ii) = t1b;
    sigma_mu2_all(ii) = t2b;

    % Triggering 
    t1g = (2*xs(5*p_tr+q_tr+1))^(-0.5);
    t2g = (2*xs(5*p_tr+2*q_tr+1))^(-0.5);

    theta_all(ii)   = xs(5*p_tr+1)*2*pi*t1g*t2g / xs(5*p_tr+3*q_tr+1);
    omega_all(ii)   = xs(5*p_tr+3*q_tr+1);
    sigma_x_all(ii) = t1g;
    sigma_y_all(ii) = t2g;
end

% Means & variances of transformed estimates
mu_bar1_mean = mean(mu_bar1_all);   
mu_bar1_var = var(mu_bar1_all,0);
mu_bar2_mean = mean(mu_bar2_all);   
mu_bar2_var = var(mu_bar2_all,0);
c_mean       = mean(c_all);         
c_var       = var(c_all,0);
d_mean       = mean(d_all);         
d_var       = var(d_all,0);
sigma_mu1_mean = mean(sigma_mu1_all); 
sigma_mu1_var = var(sigma_mu1_all,0);
sigma_mu2_mean = mean(sigma_mu2_all); 
sigma_mu2_var = var(sigma_mu2_all,0);

theta_mean   = mean(theta_all);     
theta_var   = var(theta_all,0);
omega_mean   = mean(omega_all);     
omega_var   = var(omega_all,0);
sigma_x_mean = mean(sigma_x_all);   
sigma_x_var = var(sigma_x_all,0);
sigma_y_mean = mean(sigma_y_all);   
sigma_y_var = var(sigma_y_all,0);


fprintf('\n=========== Transformed estimates across %d simulations ===========\n', n_kept);
fprintf('Background:\n');
fprintf('  mu_bar1: mean = %.6f, var = %.6e\n', mu_bar1_mean, mu_bar1_var);
fprintf('  mu_bar2: mean = %.6f, var = %.6e\n', mu_bar2_mean, mu_bar2_var);
fprintf('  c      : mean = %.6f, var = %.6e\n', c_mean,       c_var);
fprintf('  d      : mean = %.6f, var = %.6e\n', d_mean,       d_var);
fprintf('  sigma_mu1: mean = %.6f, var = %.6e\n', sigma_mu1_mean, sigma_mu1_var);
fprintf('  sigma_mu2: mean = %.6f, var = %.6e\n', sigma_mu2_mean, sigma_mu2_var);

fprintf('Triggering:\n');
fprintf('  theta  : mean = %.6f, var = %.6e\n', theta_mean,   theta_var);
fprintf('  omega  : mean = %.6f, var = %.6e\n', omega_mean,   omega_var);
fprintf('  sigma_x: mean = %.6f, var = %.6e\n', sigma_x_mean, sigma_x_var);
fprintf('  sigma_y: mean = %.6f, var = %.6e\n', sigma_y_mean, sigma_y_var);
fprintf('===================================================================\n');



fprintf('\n=========== Transformed estimates across %d simulations ===========\n', n_kept);

