clearvars
clc
%rng 4567
BIC_5 = load('BICBURG2002.mat');

params = BIC_5.BIC{1,2};

% extrating parameter for mu
g_par = reshape(params(6:13),[4,2]);

% extrating parameter for mu
mu_par = params(1:5); % q=1



mu_b0 = mu_par(1);
mu_b1 = mu_par(2);
mu_b2 = mu_par(3);
mu_b3 = mu_par(4);
mu_b4 = mu_par(5);

g_b0 = g_par(1,:);
g_b1 = g_par(2,:);
g_b2 = g_par(3,:);
g_b3 = g_par(4,:);



% disp(g_params);
mu0 = @(b0,b1,b2,b3,b4,u,v)b0*exp(-b2*(u-b1).^2-b4*(v-b3).^2);
mu_cor = -2:0.01:40;
mu1_cor = -8:0.01:40;
nor_mu = 0;

mu_x = mu0(mu_b0,mu_b1,mu_b2,mu_b3,mu_b4,mu_cor,mu_b3)*(pi/mu_b4)^0.5;
mu_y = mu0(mu_b0,mu_b1,mu_b2,mu_b3,mu_b4,mu_b1,mu1_cor)*(pi/mu_b2)^0.5;
nor_mu = nor_mu + pi*mu_b0/(mu_b2*mu_b4)^0.5;

nor_g = 0;
g0 = @(b0,b1,b2,b3,u,v,s)b0*exp(-b1*u.^2-b2*v.^2-b3*s);
g_cor = 0:0.1:20;
g1_cor = -1:0.02:1;
g2_cor = -4:0.1:4;
n = length(g_cor);
n1 = length(g1_cor);
n2 = length(g2_cor);
g_x = zeros(2,n1);
g_y = zeros(2,n2);
g_t = zeros(2,n);
for i=1:2
g_x(i,:) = g0(g_b0(i),g_b1(i),g_b2(i),g_b3(i),g1_cor,0,0)*(pi/g_b2(i))^0.5/g_b3(i);
g_y(i,:) = g0(g_b0(i),g_b1(i),g_b2(i),g_b3(i),0,g2_cor,0)*(pi/g_b1(i))^0.5/g_b3(i);
g_t(i,:) = g0(g_b0(i),g_b1(i),g_b2(i),g_b3(i),0,0,g_cor)*pi/(g_b1(i)*g_b2(i))^0.5;
nor_g =  nor_g + pi*g_b0(i)/(g_b1(i)*g_b2(i))^0.5/g_b3(i);
end


figure(1)
subplot(2,3,1)
plot(mu_cor,sum(mu_x,1)/nor_mu)
xticks(-2:10:40)

subplot(2,3,2)
plot(mu1_cor,sum(mu_y,1)/nor_mu)
xticks(-20:10:40)
%subplot(2,3,3)
%plot(mu1_cor,sum(mu_t,1)/nor_mu)
subplot(2,3,4)
semilogy(g1_cor,sum(g_x,1)/nor_g)
hold on
ylims = ylim; %[min(sum(g_x,1)/nor_g),max(sum(g_x,1)/nor_g)];
line([-0.55, -0.55],ylims,'color','r', 'Linestyle','--');
line([0.55, 0.55],ylims,'color','r', 'Linestyle','--');
text(-0.55,mean(ylims),'-0.55','Fontsize',10,'HorizontalAlignment', 'right');
text(0.55,mean(ylims),'0.55','Fontsize',10);
xticks(-1:0.4:1);
xlabel('Km')
hold off
%xlim([-6,6]);
%semilogy(x,y);
subplot(2,3,5)
semilogy(g2_cor,sum(g_y,1)/nor_g)
hold on
ylim = [min(sum(g_y,1)/nor_g),max(sum(g_y,1)/nor_g)];
line([-1.8, -1.8],ylim,'color','r', 'Linestyle','--');
line([1.8, 1.8],ylim,'color','r', 'Linestyle','--');
text(-1.8,mean(ylim),'-1.8','Fontsize',10,'HorizontalAlignment', 'right');
text(1.8,mean(ylim),'1.8','Fontsize',10);
xticks(-5:1:5);
xlabel('Km')
hold off
%xlim([-6,6]);
%semilogy(x,y);
subplot(2,3,6)
plot(g_cor,sum(g_t,1)/nor_g)
xticks(0:4:20);

xlabel('Time(days)')
