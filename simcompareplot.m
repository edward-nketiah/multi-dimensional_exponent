clearvars, clc

% the background function
mu0 = @(b0,b1,b2,b3,b4,u,v)b0*exp(-b2*(u-b1).^2-b4*(v-b3).^2);
mu_cor = -10:0.5:30;
mu1_cor = 0.1:0.1:100;

% triggering function
g0 = @(b0,b1,b2,b3,u,v,s)b0*exp(-b1*u.^2-b2*v.^2-b3*s);
g_cor = 0:0.5:40;
g1_cor = -0.05:0.001:0.05;
g2_cor = -0.5:0.01:0.5;


% estimated values for the 5 x 10^3 simulated data
mu5_b0 =  0.0463844638673967;%77764200233;a_j
mu5_b1 = 9.89566447989127;  % c
mu5_b2 = 0.0251560510746204;%91358024691; alpha_j
mu5_b3 = 10.1284501946275;  % d
mu5_b4 = 0.0254414570887354;%91358024691; beta_j

g5_b0 = 3.92973261510817;%98861837907; b_j
g5_b1 = 5100.24377451116; % alpha_j
g5_b2 = 51.7494960110381;   % beta_j
g5_b3 = 0.125798286764632;  % omega

% these functions are for the 5x10^3 values plot
mu_x = mu0(mu5_b0,mu5_b1,mu5_b2,mu5_b3,mu5_b4,mu_cor,mu5_b3)*(pi/mu5_b4)^0.5;
mu_y = mu0(mu5_b0,mu5_b1,mu5_b2,mu5_b3,mu5_b4,mu5_b1,mu_cor)*(pi/mu5_b2)^0.5;
%mu_t = mu0(mu5_b0,mu5_b1,mu5_b2,mu5_b3,mu5_b4,mu5_b1,0,mu1_cor)*pi/(mu5_b2*mu5_b4)^0.5;
nor_mu = pi*mu5_b0/(mu5_b2*mu5_b4)^0.5;

g_x = g0(g5_b0,g5_b1,g5_b2,g5_b3,g1_cor,0,0)*(pi/g5_b2)^0.5/g5_b3;
g_y = g0(g5_b0,g5_b1,g5_b2,g5_b3,0,g2_cor,0)*(pi/g5_b1)^0.5/g5_b3;
g_t = g0(g5_b0,g5_b1,g5_b2,g5_b3,0,0,g_cor)*pi/(g5_b1*g5_b2)^0.5;
nor_g =  pi*g5_b0/(g5_b1*g5_b2)^0.5/g5_b3;


% estimated values for the 5 x 10^4 simulated data
mu10_b0 = 0.0450310087850960;%77764200233;a_j
mu10_b1 = 10.0383286501518;  % c
mu10_b2 = 0.0243472701410085;%91358024691; alpha_j
mu10_b3 = 10.0214531721950;  % d
mu10_b4 = 0.0252763492505892;%91358024691; beta_j

g10_b0 = 3.18681841799280;%98861837907; b_j
g10_b1 = 5122.35088994370; % alpha_j
g10_b2 = 50.0081694767523;   % beta_j
g10_b3 = 0.0998251216799322;  % omega

% these functions are for the 1 x10^5 values plot
mu_x2 = mu0(mu10_b0,mu10_b1,mu10_b2,mu10_b3,mu10_b4,mu_cor,mu10_b3)*(pi/mu10_b4)^0.5;
mu_y2 = mu0(mu10_b0,mu10_b1,mu10_b2,mu10_b3,mu10_b4,mu10_b1,mu_cor)*(pi/mu10_b2)^0.5;
%mu_t2 = mu0(mu10_b0,mu10_b1,mu10_b2,mu10_b3,mu10_b4,mu10_b5,mu10_b1,0,mu1_cor)*pi/(mu10_b2*mu10_b4)^0.5;
nor_mu2 = pi*mu10_b0/(mu10_b2*mu10_b4)^0.5;

g_x2 = g0(g10_b0,g10_b1,g10_b2,g10_b3,g1_cor,0,0)*(pi/g10_b2)^0.5/g10_b3;
g_y2 = g0(g10_b0,g10_b1,g10_b2,g10_b3,0,g2_cor,0)*(pi/g10_b1)^0.5/g10_b3;
g_t2 = g0(g10_b0,g10_b1,g10_b2,g10_b3,0,0,g_cor)*pi/(g10_b1*g10_b2)^0.5;
nor_g2 =  pi*g10_b0/(g10_b1*g10_b2)^0.5/g10_b3;


% these are values obtain after converting the true values to the 
% same range as the estimated using calculator2.m file 
mu_b0 = 0.0449;%
mu_b1 = 10;
mu_b3 = 10;
mu_b2 = 0.0247;%
mu_b4 = 0.0247;%
mu_b5 = 0.1;
g_b0 = 3.1831;%
g_b1 = 5000;
g_b2 = 50;
g_b3 = 0.1;

% these functions are for the true values plot
mu_x1 = mu0(mu_b0,mu_b1,mu_b2,mu_b3,mu_b4,mu_cor,mu_b3)*(pi/mu_b4)^0.5;
mu_y1 = mu0(mu_b0,mu_b1,mu_b2,mu_b3,mu_b4,mu_b1,mu_cor)*(pi/mu_b2)^0.5;
%mu_t1 = mu0(mu_b0,mu_b1,mu_b2,mu_b3,mu_b4,mu_b5,mu_b1,0,mu1_cor)*pi/(mu_b2*mu_b4)^0.5;
nor_mu1 = pi*mu_b0/(mu_b2*mu_b4)^0.5;

g_x1 = g0(g_b0,g_b1,g_b2,g_b3,g1_cor,0,0)*(pi/g_b2)^0.5/g_b3;
g_y1 = g0(g_b0,g_b1,g_b2,g_b3,0,g2_cor,0)*(pi/g_b1)^0.5/g_b3;
g_t1 = g0(g_b0,g_b1,g_b2,g_b3,0,0,g_cor)*pi/(g_b1*g_b2)^0.5;
nor_g1 =  pi*g_b0/(g_b1*g_b2)^0.5/g_b3;


figure(1)
subplot(2,3,1)
plot(mu_cor,sum(mu_x,1)/nor_mu,'o')
hold on
plot(mu_cor,sum(mu_x2,1)/nor_mu2,'k+')
hold on
plot(mu_cor,sum(mu_x1,1)/nor_mu1,'r-')
xlabel('\mu(x)')
subplot(2,3,2)
plot(mu_cor,sum(mu_y,1)/nor_mu,'o')
hold on 
plot(mu_cor,sum(mu_y2,1)/nor_mu2,'k+')
hold on
plot(mu_cor,sum(mu_y1,1)/nor_mu1,'r-')
xlabel('\mu(y)')
%subplot(2,3,3)
%plot(mu1_cor,sum(mu_t,1)/nor_mu,'bo')
%hold on 
%plot(mu1_cor,sum(mu_t1,1)/nor_mu1)
subplot(2,3,4)
plot(g1_cor,sum(g_x,1)/nor_g,'o')
hold on 
plot(g1_cor,sum(g_x2,1)/nor_g2,'k+')
hold on 
plot(g1_cor,sum(g_x1,1)/nor_g1,'r-')
xlabel('g(x)')
subplot(2,3,5)
plot(g2_cor,sum(g_y,1)/nor_g,'o')
hold on 
plot(g2_cor,sum(g_y2,1)/nor_g2,'k+')
hold on 
plot(g2_cor,sum(g_y1,1)/nor_g1,'r-')
xlabel('g(y)')
subplot(2,3,6)
plot(g_cor,sum(g_t,1)/nor_g,'o')
hold on
plot(g_cor,sum(g_t2,1)/nor_g2,'k+') 
hold on
plot(g_cor,sum(g_t1,1)/nor_g1,'r-')
xlabel('g(t)')
legend('Estimated values with 5 x 10^3 data size','Estimated values with 5 x 10^4 data size','True Values')
