clear
clc
L5 = load('L_51.mat');
L5 = L5.L;
L1_5 = load('L1_51.mat');
L1_5 = L1_5.L1;
dist5 = load('dist_51');
dist5 = dist5.dist;


L10 = load('L_50k.mat');
L10 = L10.L;
L1_10 = load('L1_50k.mat');
L1_10 = L1_10.L1;
dist10 = load('dist_50k.mat');
dist10 = dist10.dist;
figure(1)
subplot(1, 2, 1)
plot(dist5,min(L5),'b-')
hold on 
plot(dist5,max(L5),'b-')
hold on
shadedplot(dist5,min(L5),max(L5),'c');
grid on
hold on
plot(dist5,mean(L5),'k-')
%plot(dist,upper_bound_L,'r-')
%hold on
plot(dist5,[min(L1_5);max(L1_5)],'r-')
hold off
xlabel('d')
ylabel('L(d)-d')

subplot(1, 2, 2)
plot(dist10,min(L10),'b-')
hold on 
plot(dist10,max(L10),'b-')
hold on
shadedplot(dist10,min(L10),max(L10),'c');
grid on
hold on
plot(dist10,mean(L10),'k-')
%plot(dist,upper_bound_L,'r-')
%hold on
plot(dist10,[min(L1_10);max(L1_10)],'r-')
hold off
xlabel('d')
ylabel('L(d)-d')