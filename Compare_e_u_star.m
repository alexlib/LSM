%following Shao 2013
clear all; close all; clc;
% formatting
addpath('/Users/travismorrison/Documents/Code/functions_library')
ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);
% looking at Playa on May 25th
load('/Users/travismorrison/Documents/Local_Data/MATERHORN/data/tower_data/Playa_tower_raw/PlayaSpring_raw_GPF_LinDet_2013_05_24.mat')
start_index = 432000; %May 24 0600 UTC or 0000MST
end_index = start_index +(20*60*60*24); %May 25 0600 UTC or 0200MST
% load radiation data: 5 min data
data_start = 6695; %24 May 0600 UTC
data_end = data_start + (12*24); %25 May 0600 UTC
load('./Materhorn_data/MATERHORN_Rad_data.mat');
SWdn = interp1(linspace(0,24,289),rad_data(data_start:data_end,6),linspace(0,24,1440),'spline');
LWup = interp1(linspace(0,24,289),rad_data(data_start:data_end,9),linspace(0,24,1440),'spline');
clear rad_data
%
u = rearrangeHeights(rawFlux.uPF(start_index:end_index,:));
v = rearrangeHeights(rawFlux.vPF(start_index:end_index,:));
w = rearrangeHeights(rawFlux.wPF(start_index:end_index,:));
theta = rearrangeHeights(rawFlux.fwTh(start_index:end_index,:));


u_prime = rearrangeHeights(rawFlux.uPF_Prime(start_index:end_index-1,:));
v_prime = rearrangeHeights(rawFlux.vPF_Prime(start_index:end_index-1,:));
w_prime = rearrangeHeights(rawFlux.wPF_Prime(start_index:end_index-1,:));
theta_prime = rearrangeHeights(rawFlux.fwThPrime(start_index:end_index-1,:));
clear rawFlux


chunk = 20*60*10;
%compute terms on 15 min chunk
u_prime_sq = squeeze(mean(reshape((u_prime.^2)',[chunk,144,6]),1,'omitnan'));
v_prime_sq = squeeze(mean(reshape((v_prime.^2)',[chunk,144,6]),1,'omitnan'));
w_prime_sq = squeeze(mean(reshape((w_prime.^2)',[chunk,144,6]),1,'omitnan'));

%% compute turbulent Pr #

wTheta_prime = w_prime.*theta_prime;

%compute derivatives
dt = 1/20;
du_dt = (u(2:end,:)-u(1:end-1,:))./dt;

dz_2 = 2.01-1.99;
P = zeros(501,1);
dT_dz = zeros(1728001,1);
for i = 1:(end_index-start_index)
    %fit to poly
    P(:) = lagrangepoly([5,2.02,0.061],theta(i,4:6),0:.01:5);
    %take derivate as a central difference around 2.02 m. There for with
    %interpalation of 501 points we have a 2dz = .02
    %
    dT_dz(i) = (P(204) - P(202))/dz_2;
    if mod(i,10000)==0
        index = i
    end
end


% Calculate the unfiltered variables 
uw_prime = squeeze(mean(reshape((u_prime.*w_prime)',[chunk,144,6]),1,'omitnan'));
vw_prime = squeeze(mean(reshape((v_prime.*w_prime)',[chunk,144,6]),1,'omitnan'));

%Compute the turbulent Pr number then Chunk it %%WRONG%%, should be du_dz
Pr = squeeze(mean(reshape(((u_prime(1:end,5).*w_prime(1:end,5).*dT_dz(1:end-1))./...
    (wTheta_prime(1:end,5).*du_dt(1:end,5))),[chunk,144]),1,'omitnan'));
%Compute tke and u_star
tke = 0.5.*(u_prime_sq+v_prime_sq+...
    w_prime_sq);

u_star = (uw_prime.^2+vw_prime.^2).^(1/4);

% plot tke and u_star
figure()
plot(linspace(0,24,144),tke(:,5).^.5,'k-*')
hold on
plot(linspace(0,24,144),u_star(:,5),'k-o')
legend('$\sqrt{e}$ [m$^2$s$^{-2}$]','$u^*$ [ms$^{-1}$]')
xlabel("time [hrs]")
axis tight
grid on
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);
% Place second set of axes on same plot
handaxes2 = axes('Position', [0.7 0.5 0.2 0.3]);

plot(linspace(0,24,1440),SWdn,  'k--')
%set(handaxes2, 'Box', 'off')
xlabel('t (hrs)')
ylabel('SWdn [Wm$^{-2}$]')
axis tight
grid on

%% Use taylors frozen theory hypothesis
load('./Materhorn_data/playaSpring30minLinDetUTESpac3.mat')
%load tower data: 30 min data
start_index = 1068;%24 May 0600 UTC
end_index = start_index + (2*24);%25 May 0600 UTC
U = rearrangeHeights(playaSpring.spdAndDir(start_index:end_index, [3:3:18]));
L = rearrangeHeights(playaSpring.L(start_index:end_index,2:end));
%clear playaSpring

freq = 20;
npoints = size(U,1);
chunk = freq*60*30; %in frames
time = chunk/freq; % [sec]
tower_height = 5;
les_grid_size = 10; %[m]

t = (U(:,tower_height)./les_grid_size).^(-1);

%
% figure()
% plot(linspace(0,24,49),t,'k-')
% hold on
% ylabel('$\frac{U}{dx_{LES}} $')
% xlabel('time [hrs]')

filter_size = t.*freq; %[npoints for flux data]
taylors_window = chunk; %[npoints = to length of application of taylors hypothesis]


%% Filter data into Resolved and SGS
tmp =  1728000-sum(floor(filter_size))+8;
u_prime_resolved = zeros(1,tmp);
v_prime_resolved = zeros(1,tmp);
w_prime_resolved = zeros(1,tmp);
theta_prime_resolved= zeros(1,tmp);
u_resolved = zeros(tmp,6);


start_index = 1;
nchunks = 48;
for i = 1:nchunks
    end_index = floor(start_index + filter_size(i));
    for j = 1:taylors_window-filter_size(i)
        
        %since we need theta for the finite difference in height
        
        %primes
        u_prime_resolved(start_index) = squeeze(mean(u_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        v_prime_resolved(start_index) = squeeze(mean(v_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        w_prime_resolved(start_index) = squeeze(mean(w_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        theta_prime_resolved(start_index) = squeeze(mean(theta_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        start_index = start_index+1;
        end_index = end_index+1;
    end
    if mod(i,1)==0
          chunk_num = i
    end
end
%
%% perform the filtering at each height for theta since we have to take the
% vertical derivative
start_index = 1;
nchunks = 48;
theta_resolved = zeros(tmp,6);
for k = 1:6
    start_index = 1;
    for i = 1:nchunks
        end_index = floor(start_index + filter_size(i));
        for j = 1:taylors_window-filter_size(i)
            theta_resolved(start_index,k) = squeeze(mean(theta(start_index:end_index,k),1,'omitnan')); %this solves for the grid resolved
            u_resolved(start_index,k) = squeeze(mean(u(start_index:end_index,k),1,'omitnan')); %this solves for the grid resolved
            start_index = start_index+1;
            end_index = end_index+1;
        end
         if mod(i,1)==0
          chunk_num = i
        end
    end
end
%% Re distribute points over which 
u_prime_resolved_interp = interp1(linspace(0,24,tmp),u_prime_resolved,linspace(0,24,1728000),'spline');
v_prime_resolved_interp = interp1(linspace(0,24,tmp),v_prime_resolved,linspace(0,24,1728000),'spline');
w_prime_resolved_interp = interp1(linspace(0,24,tmp),w_prime_resolved,linspace(0,24,1728000),'spline');
theta_prime_resolved_interp = interp1(linspace(0,24,tmp),theta_prime_resolved,linspace(0,24,1728000),'spline');
%
for k = 1:6
    theta_resolved_interp(:,k) = interp1(linspace(0,24,tmp),theta_resolved(:,k),linspace(0,24,1728000),'spline');
    u_resolved_interp(:,k) = interp1(linspace(0,24,tmp),u_resolved(:,k),linspace(0,24,1728001),'spline');
end
%
u_prime_sgs = u_prime(:,tower_height)' - u_prime_resolved_interp;
v_prime_sgs = v_prime(:,tower_height)' - v_prime_resolved_interp;
w_prime_sgs = w_prime(:,tower_height)' - w_prime_resolved_interp;
theta_prime_sgs = theta_prime(:,tower_height)' - theta_prime_resolved_interp;

%% need theta at each height
for k = 1:6
    theta_sgs(:,k) = theta(1:end-1,k) - theta_resolved_interp(:,k);
    u_sgs(:,k) = u(1:end,k) - u_resolved_interp(:,k);
end

%%

uw_prime_sgs = u_prime_sgs.*w_prime_sgs;
wTheta_prime_sgs = w_prime_sgs.*theta_prime_sgs;

%compute derivatives
dt = 1/20;
du_dt = (u_sgs(2:end,tower_height)-u_sgs(1:end-1,tower_height))./dt;
%
dz_2 = 2.01-1.99;
P = zeros(501,1);
dT_dz = zeros(1728001,1);
for i = 1:length(theta_sgs(:,1))
    %fit to poly
    P(:) = lagrangepoly([5,2.02,0.061],theta_sgs(i,4:6),0:.01:5);
    %take derivate as a central difference around 2.02 m. There for with
    %interpalation of 501 points we have a 2dz = .02
    %
    dT_dz(i) = (P(204) - P(202))/dz_2;
    if mod(i,10000)==0
        index = i
    end
end

%% Compute du dz at 2m
dz_2 = 2.01-1.99;
P = zeros(501,1);
dU_dz = zeros(1728001,1);
for i = 1:length(theta_sgs(:,1))
    %fit to poly
    P(:) = lagrangepoly([5,2.02,0.061],u_sgs(i,4:6),0:.01:5);
    %take derivate as a central difference around 2.02 m. There for with
    %interpalation of 501 points we have a 2dz = .02
    %
    dU_dz(i) = (P(204) - P(202))/dz_2;
    if mod(i,10000)==0
        index = i
    end
end

%% Compute the turbulent Pr number then Chunk it
chunk = 20*60*10;






Pr_sgs = (squeeze(mean(reshape(uw_prime_sgs,[chunk,144]),1,'omitnan')).*...
    squeeze(mean(reshape(dT_dz(1:end-1)',[chunk,144]),1,'omitnan')))./...
    (squeeze(mean(reshape(wTheta_prime_sgs,[chunk,144]),1,'omitnan')).*...
    squeeze(mean(reshape(dU_dz(1:end-1)',[chunk,144]),1,'omitnan')));

%% look at each term of the Pr sgs to see why it goes negative
figure()
subplot(2,1,1)
plot(linspace(0,24,144),squeeze(mean(reshape((uw_prime_sgs),[chunk,144]),1,'omitnan')));
hold on
plot(linspace(0,24,144),squeeze(mean(reshape(dT_dz(1:end-1)',[chunk,144]),1,'omitnan')));
plot(linspace(0,24,144),squeeze(mean(reshape(wTheta_prime_sgs,[chunk,144]),1,'omitnan')));
plot(linspace(0,24,144),squeeze(mean(reshape(dU_dz(1:end-1),[chunk,144]),1,'omitnan')));
%plot(linspace(0,24,144),squeeze(mean(reshape(w_prime_sgs,[chunk,144]),1,'omitnan')));
legend('$u^\prime w^\prime$','$\partial \theta / \partial z$',...
    '$w^\prime \theta^\prime$','$\partial u / \partial z$','$w^\prime$')
grid on
%plot(linspace(0,24,2),[0 0],'k-')
axis tight
xlabel('time (hrs)')

subplot(2,1,2)
plot(linspace(0,24,144),squeeze(mean(reshape((-uw_prime_sgs./dU_dz(1:end-1)'),[chunk,144]),1,'omitnan')));
hold on
plot(linspace(0,24,144),squeeze(mean(reshape(-wTheta_prime_sgs./dT_dz(1:end-1)',[chunk,144]),1,'omitnan')));
plot(linspace(0,24,144),Pr_sgs)
legend('$k_m = \frac{-u^\prime w^\prime}{\partial u / \partial z}$',...
    '$k_h = \frac{-w^\prime T^\prime}{\partial \theta / \partial z}$',...
    '$Pr = \frac{k_m}{k_h}$')
grid on
%plot(linspace(0,24,2),[0 0],'k-')
axis tight
xlabel('time (hrs)')


%%

g = 9.8; 
%Compute Rf #
Rf = squeeze(mean(reshape((g.*wTheta_prime_sgs)./(theta_sgs(:,5)'.*uw_prime_sgs.*dU_dz(1:end-1)'),[chunk,144]),1,'omitnan'));

%Compute tke and u_star


%compute terms on 30 min chunk
u_prime_sq_sgs = squeeze(mean(reshape((u_prime_sgs.^2)',[chunk,144]),1,'omitnan'));
v_prime_sq_sgs = squeeze(mean(reshape((v_prime_sgs.^2)',[chunk,144]),1,'omitnan'));
w_prime_sq_sgs = squeeze(mean(reshape((w_prime_sgs.^2)',[chunk,144]),1,'omitnan'));

tke_sgs = 0.5.*(u_prime_sq_sgs+v_prime_sq_sgs+...
    w_prime_sq_sgs);
%ustar vars
uw_sgs = squeeze(mean(reshape((u_prime_sgs.*w_prime_sgs)',[chunk,144]),1,'omitnan'));
vw_sgs = squeeze(mean(reshape((v_prime_sgs.*w_prime_sgs)',[chunk,144]),1,'omitnan'));

u_star_sgs = (uw_sgs.^2+vw_sgs.^2).^(1/4);

%% plot e and u star sgs and Pr
figure()
plot(linspace(0,24,144),tke_sgs(:).^(.5),'k-*')
hold on
plot(linspace(0,24,144),u_star_sgs(:),'k-o')
plot(linspace(0,24,144),Pr_sgs,'k--')
legend('$\sqrt{e}_{sgs}$ [m$^2$s$^{-2}$]','$u^*_{sgs}$ [ms$^{-1}$]','$Pr_{sgs}$')
xlabel("time [hrs]")
axis tight
ylim([-10 10])
grid on
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);
% Place second set of axes on same plot
% handaxes2 = axes('Position', [0.2 0.5 0.2 0.3]);
% 
% plot(linspace(0,24,1440),SWdn,  'k--')
% %set(handaxes2, 'Box', 'off')
% xlabel('t (hrs)')
% ylabel('SWdn [Wm$^{-2}$]')
% axis tight
% grid on
%% plot Pr vs z/L
L_interp = interp1(linspace(0,24,49),L(:,tower_height),linspace(0,24,144));
figure()
loglog((2.02./L_interp),Pr_sgs,'ko')
% hold on 
% loglog(10.^[0 10e5], [1 1],'k-')
% loglog([1 1],10.^[0 10e5], 'k-')
xlim([10e-3 10e-1])
ylim([10e-3 10e1])
ylabel('$Pr_{sgs}$')
xlabel('$\zeta = z/L$')
grid on
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);

%% plot Pr vs Rf
figure()
loglog(Rf,Pr_sgs,'ko')
% hold on 
% plot([-10e5 10e5], [0 0],'k-')
% plot([0 0],[-10e5 10e5], 'k-')
% xlim([-10 10])
% ylim([-10 10])
ylabel('$Pr_{sgs}$')
xlabel('$R_f$')
grid on
axis tight
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);
%% plot Shao and MOST vs Rf
Ck = 0.15;
kappa = 0.4;
figure()
loglog(Rf,((Ck.*tke_sgs.^(.5))./Pr_sgs),'ko')
hold on 
loglog(Rf,(kappa.*u_star_sgs),'k*')
% plot([-10e5 10e5], [0 0],'k-')
% plot([0 0],[-10e5 10e5], 'k-')
% xlim([-10 10])
% ylim([-10 10])
%ylabel('$Pr_{sgs}$')
legend('$C_k \sqrt{e_{sgs}} Pr^{-1}$','$\kappa u^*_{sgs}$')
xlabel('$R_f$')
grid on
axis tight
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);
%% plot Shao and MOST vs zeta
figure()
semilogy((2.02./L_interp),((Ck.*tke_sgs.^(.5))./Pr_sgs),'ko')
hold on 
semilogy((2.02./L_interp),(kappa.*u_star_sgs),'k*')
% plot([-10e5 10e5], [0 0],'k-')
% plot([0 0],[-10e5 10e5], 'k-')
% xlim([-10 10])
% ylim([-10 10])
%ylabel('$Pr_{sgs}$')
legend('$C_k \sqrt{e_{sgs}} Pr^{-1}$','$\kappa u^*_{sgs}$')
xlabel('$\zeta = z/L$')
grid on
axis tight
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);



%% Plot kappa ustar against tke/Pr
kappa =0.4;
Ck= 0.15;
min = 0; 
max = .4;
figure()
plot(((Ck.*tke_sgs.^(.5))./1),(kappa.*u_star_sgs),'ko')
hold on
plot(linspace(-10000,10000,2),linspace(-10000,10000,2),'k--')
title('Pr  = 1','interpreter','latex')
ylim([min  max])
xlim([min max])
xlabel("")
%axis tight
grid on

ylabel('$\kappa u^*_{sgs}$')
xlabel('$C_k \sqrt{e_{sgs}} Pr^{-1}$')
%% plot coeffs vs time
figure()
plot(linspace(0,24,144),((Ck.*tke_sgs.^(.5))./Pr_sgs),'k-o')
hold on 
plot(linspace(0,24,144),(kappa.*u_star_sgs),'k-*')
grid on 
axis tight
legend('$C_k \sqrt{e_{sgs}} Pr^{-1}$','$\kappa u^*_{sgs}$')
xlabel('time (hrs)')
%% Compute pseudo r_a with tower data
emiss = 0.93;
sigma = 5.67e-8;
z = 2.02;
z0 = 0.11e-3; %roughness length [m] for Playa (Jensen et al)
start_index = 1068;%24 May 0600 UTC
end_index = start_index + (2*24);%25 May 0600 UTC
T_air_tower = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.derivedT(start_index:end_index,2:4:22)),linspace(0,24,144),'spline');
Ts = interp1(linspace(0,24,1440),(LWup./(emiss*sigma)).^(1/4),linspace(0,24,144),'spline')-273.15;
zeta = z./L_interp(:);

ra_shao = ra_fun(z,tke_sgs,z0,'Shao',u_star_sgs,zeta',Pr_sgs);
ra_most = ra_fun(z,tke_sgs,z0,'MOST',u_star_sgs,zeta',Pr_sgs);
ra_truth =(squeeze(mean(reshape(wTheta_prime_sgs,[chunk,144]),1,'omitnan'))'./(Ts'-T_air_tower(:,5))).^(-1);

%% Plot ra vs time
figure()
plot(linspace(0,24,144),ra_shao,'k-*')
hold on
plot(linspace(0,24,144),ra_most,'k-o')
plot(linspace(0,24,144),ra_truth,'k-^')
legend('$Shao$','$MOST$','$Truth$')
ylabel('$ra_{sgs}$')
xlabel('time (hrs)')
grid on
%%

figure()
plot(linspace(0,24,144),zeta,'k-*')
hold on
plot([-0 24],[0 0],'k-')
ylabel('$\zeta = z/L$')
xlabel('time (hrs)')
grid on
axis tight
%%


min =10e-1;
max = 10e4;
figure()
subplot(3,1,1)
loglog(ra_shao,ra_most,'k*')
hold on
loglog(10.^[-1 5],10.^[-1 5],'k--')
grid on
xlim([min max])
ylim([min max])
xlabel('$ra_{Shao}$')
ylabel('$ra_{MOST}$')

subplot(3,1,2)
loglog(ra_truth,ra_shao,'k*')
grid on
hold on
loglog(10.^[0 5],10.^[0 5],'k--')
% plot(linspace(-10000,10000,2),linspace(-10000,10000,2),'k--')
% plot([0 0],linspace(-10000,10000,2),'k-')
% plot(linspace(-10000,10000,2),[0 0],'k-')
grid on
xlim([min max])
ylim([min max])
xlabel('$ra_{tower}$')
ylabel('$ra_{Shao}$')

subplot(3,1,3)
loglog(ra_truth,ra_most,'k*')
grid on
hold on
loglog(10.^[0 5],10.^[0 5],'k--')
% plot(linspace(-10000,10000,2),linspace(-10000,10000,2),'k--')
% plot([0 0],linspace(-10000,10000,2),'k-')
% plot(linspace(-10000,10000,2),[0 0],'k-')
grid on
xlim([min max])
ylim([min max])
xlabel('$ra_{tower}$')
ylabel('$ra_{MOST}$')

