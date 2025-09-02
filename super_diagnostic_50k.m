clear
clc
% Generate artificial data
tic
dbstop if error
S = [20,20];
T = [0,1350];%At 1350, there are about 5,000 samples; at 15,500, there are about 100,000 samples
[space_time_data,N1,N_b] = fgenerate_data_CRIME(S,T(2));
Bound = size(space_time_data,1);
N1 = N1(2001:Bound-2000);% we disregard the first and last 2000
N_b1 = sum(N1<=N_b);
space_time_data = space_time_data(2001:Bound-2000,:);
space_time_data(:,1:2) = space_time_data(:,1:2)+10;

x1 = space_time_data(:,1);%Space-time coordinates
y1 = space_time_data(:,2);
t1 = space_time_data(:,3);
t1 = t1-t1(1);%Time Shift
n = length(x1);
n1 = 200;%The most recent 200 points are used for calculation
X = zeros(n1-1,n-1);%This value indicates that the n1 data points closest to the time point are used to calculate g in lambda.
Y = X;
T = X;
for j = 1:n1-1
    X(j,:) = x1(2:end)-[-Inf*ones(j-1,1);x1(1:end-j)];
    Y(j,:) = y1(2:end)-[-Inf*ones(j-1,1);y1(1:end-j)];
    T(j,:) = t1(2:end)-[-Inf*ones(j-1,1);t1(1:end-j)];
end


%  BIC_5 = load('BIC5NN.mat');
% BIC_5 = BIC_5.BIC;
BIC_10 = load('BIC10NN.mat');
BIC_10 = BIC_10.BIC;

S = [0,20,0,20];
dist = 0:0.1:9;
% L1 = zeros(19,length(dist));
L = zeros(19,length(dist));
L1 = zeros(19,length(dist));
x = BIC_10{1,1};
p = 0 + 1;
q = 0 + 1;
for i = 1:19
 [data,~,k] = di_method(x,x1,y1,t1,X,Y,T,S,p,q);
 n = poissrnd(k*S(2)*S(4)*1350); 
  if n>0
    space_time_data0 = [rand(n,1)*S(2),rand(n,1)*S(4),rand(n,1)*1350];
    x10 = space_time_data0(:,1);%Space-time coordinates
    y10 = space_time_data0(:,2);
    t10 = space_time_data0(:,3);
    t10 = t10-t10(1);%Time Shift
    n0 = length(x10);
    n10 = min(200,n0);%The most recent 200 points are used for calculation
    X0 = zeros(n10-1,n0-1);%This value indicates that the n1 data points closest to the time point are used to calculate g in lambda.
    Y0 = X0;
    T0 = X0;
    for j = 1:n10-1
      X0(j,:) = x10(2:end)-[-Inf*ones(j-1,1);x10(1:end-j)];
      Y0(j,:) = y10(2:end)-[-Inf*ones(j-1,1);y10(1:end-j)];
      T0(j,:) = t10(2:end)-[-Inf*ones(j-1,1);t10(1:end-j)];
    end
    [~,lambda,~] = di_method(x,x10,y10,t10,X0,Y0,T0,S,p,q);   
    temp = min((k./lambda),1) >rand(1,n);
    % using the super thin simulate inhomogeneous Poisson process with rate
    % max(k-lambda,0)
    super_thin_rate = max((k-lambda),0);
    n_super = poissrnd(sum(super_thin_rate)*S(2)*S(4)*1350);
    sup_rand_vals = rand(n_super,3);
    super_thin_data = [sup_rand_vals(:,1)*S(2),sup_rand_vals(:,2)*S(4),sup_rand_vals(:,3)*1350];
    
    data = [space_time_data0(temp,:);super_thin_data];
  end
 K = RipleysK(data(:,1:2),dist,S);
 N = poissrnd(size(data,1));
 x0 = [rand(N,1)*20,rand(N,1)*20];
 K1 = RipleysK(x0,dist,S);
 L(i,:) = sqrt(K/pi)-dist';
 L1(i,:) = sqrt(K1/pi)-dist';
end

% Plot the L-functions and their confidence intervals
figure(1)
%subplot(1, 2, 1)

plot(dist,min(L),'b-')
hold on 
plot(dist,max(L),'b-')
hold on
shadedplot(dist,min(L),max(L),'c');
grid on
hold on
plot(dist,mean(L),'k-')
%plot(dist,upper_bound_L,'r-')
%hold on
plot(dist,[min(L1);max(L1)],'r-')
hold off


