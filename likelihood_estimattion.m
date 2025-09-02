clearvars
clc
%%Generate artificial data
tic
dbstop if error
S = [20,20];
T = [0,1350];%At 1350, there are about 5,000 samples; at 15,500, there are about 100,000 samples
[space_time_data,N1,N_b] = fgenerate_data_CRIME21(S,T(2));
% [space_time_data,N1,N_b] = fgenerate_data_ETAS(S,T(2));
Bound = size(space_time_data,1);
N1 = N1(2001:Bound-2000);
N_b1 = sum(N1<=N_b);
space_time_data = space_time_data(2001:Bound-2000,:);
space_time_data(:,1:2) = space_time_data(:,1:2)+10;
space_time_data = space_time_data(:,1:3);

x1 = space_time_data(:,1);%Space-time coordinates
y1 = space_time_data(:,2);
t1 = space_time_data(:,3);
t1 = t1-t1(1);%Time Shift
n = length(x1);
n1 = 500;%The most recent 200 points are used for calculation
X = zeros(n1-1,n-1);%This value indicates that the n1 data points closest to the time point are used to calculate g in lambda
Y = X;
T = X;
for j = 1:n1-1
    X(j,:) = x1(2:end)-[-Inf*ones(j-1,1);x1(1:end-j)];
    Y(j,:) = y1(2:end)-[-Inf*ones(j-1,1);y1(1:end-j)];
    T(j,:) = t1(2:end)-[-Inf*ones(j-1,1);t1(1:end-j)];
end

min_fval = Inf(1,9); % set the initial min fval to a large number
best_init_all = {}; % initialize to store the best_init_a
for init_a =1:10
    for i = 1:9
        p=(i>3)+(i>6);%Functions and Approximations
        q=mod(i-1,3);
        int_a = [rand(1,5*(p+1)),rand(1,4*(q+1))];
        disp(length(int_a))
        tic
        % [x, fval, exitflag, output] = fmincon(@(a)li_fun(a,space_time_data,X,Y,T,p+1,q+1),int_a,[],[],[],[],zeros(1,4*(p+q+2)),[],[]);
        [~, fval, exitflag, output] = fmincon(@(a)li_fun1(a,x1,y1,t1,X,Y,T,[0 20 0 20],p+1,q+1),int_a,[],[],[],[],zeros(1,5*(p+1)+4*(q+1)),[],@(x)mycon(x,p+1,q+1));
        toc
        %BIC{1,i} = x;
        %BIC{2,i} = fval;
        %BIC{3,i} = fval+(5*(p+1)+4*(q+1))*log(n);
        %check if current fval is min
        if fval < min_fval(i)
           min_fval(i) = fval;
           best_init_all{i} = int_a; % store the corresponding init_a
        end
    end
end
BIC=cell(3,9);
for i = 1:9
p=(i>3)+(i>6);%Functions and Approximations
q=mod(i-1,3);
best_init_a = best_init_all{i};
disp(length(best_init_a))
tic
% [x, fval, exitflag, output] = fmincon(@(a)li_fun(a,space_time_data,X,Y,T,p+1,q+1),int_a,[],[],[],[],zeros(1,4*(p+q+2)),[],[]);
[x, fval, exitflag, output] = fmincon(@(a)li_fun1(a,x1,y1,t1,X,Y,T,[0 20 0 20],p+1,q+1),best_init_a,[],[],[],[],zeros(1,5*(p+1)+4*(q+1)),[],@(x)mycon(x,p+1,q+1));
toc
BIC{1,i} = x;
BIC{2,i} = fval;
BIC{3,i} = fval+(5*(p+1)+4*(q+1))*log(n);
end 

%Plotting BIC
x = BIC{1,1};
mu0 = @(b0,b1,b2,b3,b4,u,v)b0*exp(-b2*(u-b1).^2-b4*(v-b3).^2);
nor_mu = 0;
mu_cor = -50:0.05:50;
mu1_cor = 0.1:0.1:100;
n = length(mu_cor);
n1 = length(mu1_cor);
p = 0+1;
q = 0+1;
% p = 1;
% q = 2;
mu_x = zeros(p,n);
mu_y = zeros(p,n);
%mu_t = zeros(p,n1);
for j = 1:p
  mu_x(j,:) = mu0(x(j),x(p+j),x(2*p+j),x(3*p+j),x(4*p+j),mu_cor,x(3*p+j))*(pi/x(4*p+j))^0.5;
  mu_y(j,:) = mu0(x(j),x(p+j),x(2*p+j),x(3*p+j),x(4*p+j),x(p+j),mu_cor)*(pi/x(2*p+j))^0.5;
  %mu_t(j,:) = mu0(x(j),x(p+j),x(2*p+j),x(3*p+j),x(4*p+j),x(5*p+j),x(p+j),x(3*p+j),mu1_cor)*pi/(x(2*p+j)*x(4*p+j))^0.5;
  nor_mu = nor_mu + pi*x(j)/(x(2*p+j)*x(4*p+j))^0.5;
end
nor_g = 0;
g0 = @(b0,b1,b2,b3,u,v,s)b0*exp(-b1*u.^2-b2*v.^2-b3*s);
g_cor = 0:0.01:50;
g1_cor = -0.05:0.001:0.05;
g2_cor = -0.5:0.01:0.5;

%g1_cor = -1:0.01:1;
%g2_cor = -1:0.01:1;
n = length(g_cor);
n1 = length(g1_cor);
n2 = length(g2_cor);
g_x = zeros(q,n1);
g_y = zeros(q,n2);
g_t = zeros(q,n);
for j = 1:q
    g_x(j,:) = g0(x(5*p+j),x(5*p+q+j),x(5*p+2*q+j),x(5*p+3*q+j),g1_cor,0,0)*(pi/x(5*p+2*q+j))^0.5/x(5*p+3*q+j);
    g_y(j,:) = g0(x(5*p+j),x(5*p+q+j),x(5*p+2*q+j),x(5*p+3*q+j),0,g2_cor,0)*(pi/x(5*p+q+j))^0.5/x(5*p+3*q+j);
    g_t(j,:) = g0(x(5*p+j),x(5*p+q+j),x(5*p+2*q+j),x(5*p+3*q+j),0,0,g_cor)*pi/(x(5*p+q+j)*x(5*p+2*q+j))^0.5;
    nor_g = nor_g + pi*x(5*p+j)/(x(5*p+q+j)*x(5*p+2*q+j))^0.5/x(5*p+3*q+j);
end
% these are values obtain after converting the true values to the 
% same range as the estimated using calculator2.m file 
% mu_b0 = 0.0449;%77764200233;
% mu_b1 = 10;
% mu_b3 = 10;
% mu_b2 = 0.0247;%91358024691;
% mu_b4 = 0.0247;%91358024691;
% mu_b5 = 0.1;
% g_b0 = 3.1831;%98861837907;
% g_b1 = 5000;
% g_b2 = 50;
% g_b3 = 0.1;

% these functions are for the true values plot
% mu_x1 = mu0(mu_b0,mu_b1,mu_b2,mu_b3,mu_b4,mu_cor,mu_b3)*(pi/mu_b4)^0.5;
% mu_y1 = mu0(mu_b0,mu_b1,mu_b2,mu_b3,mu_b4,mu_b1,mu_cor)*(pi/mu_b2)^0.5;
% %mu_t1 = mu0(mu_b0,mu_b1,mu_b2,mu_b3,mu_b4,mu_b5,mu_b1,0,mu1_cor)*pi/(mu_b2*mu_b4)^0.5;
% nor_mu1 = pi*mu_b0/(mu_b2*mu_b4)^0.5;
% 
% g_x1 = g0(g_b0,g_b1,g_b2,g_b3,g1_cor,0,0)*(pi/g_b2)^0.5/g_b3;
% g_y1 = g0(g_b0,g_b1,g_b2,g_b3,0,g2_cor,0)*(pi/g_b1)^0.5/g_b3;
% g_t1 = g0(g_b0,g_b1,g_b2,g_b3,0,0,g_cor)*pi/(g_b1*g_b2)^0.5;
% nor_g1 =  pi*g_b0/(g_b1*g_b2)^0.5/g_b3;

% figure(1)
% subplot(2,3,1)
% plot(mu_cor,sum(mu_x,1)/nor_mu,'bo')
% hold on 
% plot(mu_cor,sum(mu_x1,1)/nor_mu1,'r-')
% xlabel('\mu(x)')
% subplot(2,3,2)
% plot(mu_cor,sum(mu_y,1)/nor_mu,'o')
% hold on 
% plot(mu_cor,sum(mu_y1,1)/nor_mu1,'r-')
% xlabel('\mu(y)')
% subplot(2,3,3)
% plot(mu1_cor,sum(mu_t,1)/nor_mu,'bo')
% hold on 
% plot(mu1_cor,sum(mu_t1,1)/nor_mu1)
% subplot(2,3,4)
% plot(g1_cor,sum(g_x,1)/nor_g,'o')
% %hold on 
% plot(g1_cor,sum(g_x1,1)/nor_g1,'r-')
% xlabel('g(x)')
% subplot(2,3,5)
% plot(g2_cor,sum(g_y,1)/nor_g,'o')
% hold on 
% plot(g2_cor,sum(g_y1,1)/nor_g1,'r-')
% xlabel('g(y)')
% subplot(2,3,6)
% plot(g_cor,sum(g_t,1)/nor_g,'o')
% hold on
% plot(g_cor,sum(g_t1,1)/nor_g1,'r-')
% xlabel('g(t)')
% legend('Estimated values with 5 x 10^3 data size','True Values')
%'Estimated values with 1x10^5 data size')
% figure(1)
% x_data = space_time_data(space_time_data(:,4)==0,1);
% y_data = space_time_data(space_time_data(:,4)==0,2);
% plot(x_data,y_data,'.blue')
% hold on;
% x_data1 = space_time_data(space_time_data(:,4)==1,1);
% y_data1 = space_time_data(space_time_data(:,4)==1,2);
% plot(x_data1,y_data1,'.red')
% hold off;
% axis([0 20 0 20])

figure(1)
subplot(2,3,1)
plot(mu_cor,sum(mu_x,1)/nor_mu)
subplot(2,3,2)
plot(mu_cor,sum(mu_y,1)/nor_mu)
%subplot(2,3,3)
%plot(mu1_cor,sum(mu_t,1)/nor_mu)
subplot(2,3,4)
plot(g1_cor,sum(g_x,1)/nor_g)
subplot(2,3,5)
plot(g2_cor,sum(g_y,1)/nor_g)
subplot(2,3,6)
plot(g_cor,sum(g_t,1)/nor_g)