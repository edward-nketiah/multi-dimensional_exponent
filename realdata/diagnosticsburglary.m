clear
clc
%rng(100)%0110

data2 = load('Burglary2002.mat');
data3 = data2.data1;

lat = data3(:,1);

lon = data3(:,2);

% Convert latitude and longitude to kilometers
% Using a reference latitude (e.g., the mean latitude) for accurate scaling
reference_lat = mean(lat);
lat_km = (lat - reference_lat) * 111; % Latitude scale is roughly 111 km per degree
lon_km = (lon - mean(lon)) * 111 * cosd(reference_lat); % Adjust longitude for latitude

% Shift coordinates to start from zero
lat_km = lat_km - min(lat_km);
lon_km = lon_km - min(lon_km);


x1 = lon_km;
y1 = lat_km;

S = [15 23 15 30];
t1 = data3(:,3);

n = length(x1);
n1 = 200;%The most recent 200 points are used for calculation
X = zeros(n1-1,n-1);%This value indicates that the n1 data points closest to the time point are used to calculate g in lambda
Y = X;
T = X;
%eps =1e-4;
for j = 1:n1-1
    X(j,:) = x1(2:end)-[-Inf*ones(j-1,1);x1(1:end-j)];
    Y(j,:) = y1(2:end)-[-Inf*ones(j-1,1);y1(1:end-j)];
    T(j,:) = t1(2:end)-[-Inf*ones(j-1,1);t1(1:end-j)];
end


BIC_5 = load('BICBURG2002.mat'); % felony crime
BIC_5 = BIC_5.BIC;
rng(115)
dist = 0:0.01:3;
% L1 = zeros(19,length(dist));
L = zeros(3,length(dist));
L1 = zeros(3,length(dist));
x = BIC_5{1,2};
p =0+1;
q = 1+1;

for i = 1:10
 [data,~,k] = di_method(x,x1,y1,t1,X,Y,T,S,p,q);
 n = poissrnd(k*(S(2)-S(1))*(S(4)-S(3))*364.9); %length(data3)
  if n>0
    rand_vals = rand(n,3);
    space_time_data0 = [rand_vals(:,1)*(S(2)-S(1))+S(1),rand_vals(:,2)*(S(4)-S(3))+ S(3),rand_vals(:,3)*364.9];
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
    n_super = poissrnd(sum(super_thin_rate)*(S(2)-S(1))*(S(4)-S(3))*364.9);
    sup_rand_vals = rand(n_super,3);
    super_thin_data = [sup_rand_vals(:,1)*(S(2)-S(1))+S(1),sup_rand_vals(:,2)*(S(4)-S(3))+S(3),sup_rand_vals(:,3)*364.9];

    %combine the two process
    data = [space_time_data0(temp,:);super_thin_data];
  end 
    
  
 K = RipleysK(data(:,1:2),dist,S);
 N = poissrnd(size(data,1));
 x0 = [rand(N,1)*(S(2)-S(1))+S(1),rand(N,1)*(S(4)-S(3))+S(3)];
 K1 = RipleysK(x0,dist,S);
 L(i,:) = sqrt(K/pi)-dist';
 L1(i,:) = sqrt(K1/pi)-dist';
end


% Plot the L-functions and their confidence intervals
figure(1)

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



