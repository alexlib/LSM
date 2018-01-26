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


chunk = 20*60*15;
%compute terms on 15 min chunk
u_prime_sq = squeeze(mean(reshape((u_prime.^2)',[chunk,144,6]),1,'omitnan'));
v_prime_sq = squeeze(mean(reshape((v_prime.^2)',[chunk,144,6]),1,'omitnan'));
w_prime_sq = squeeze(mean(reshape((w_prime.^2)',[chunk,144,6]),1,'omitnan'));

uw_prime = squeeze(mean(reshape((u_prime.*w_prime)',[chunk,144,6]),1,'omitnan'));
vw_prime = squeeze(mean(reshape((v_prime.*w_prime)',[chunk,144,6]),1,'omitnan'));
wTheta_prime = squeeze(mean(reshape((w_prime.*theta_prime)',[chunk,144,6]),1,'omitnan'));

%compute derivatives
dt = 1/20; 
du_dt = squeeze(mean(reshape(((u(2:end,:)-u(1:end-1,:))./dt)',[chunk,144,6]),1,'omitnan')); 



%Compute tke and u_star
tke = 0.5.*(u_prime_sq+v_prime_sq+...
    w_prime_sq);

u_star = (uw_prime.^2+vw_prime.^2).^(1/4);

%% plot
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
start_index = 1072;%25 May 0000 UTC
end_index = start_index + (2*24);%26 May 0000 UTC
U = rearrangeHeights(playaSpring.spdAndDir(start_index:end_index, [3:3:18]));
clear playaSpring

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
u_resolved = zeros(1,tmp);
v_resolved = zeros(1,tmp);
w_resolved = zeros(1,tmp);
theta_resolved= zeros(1,tmp);

start_index = 1;
nchunks = 48;
for i = 1:nchunks
    end_index = floor(start_index + filter_size(i));
    for j = 1:taylors_window-filter_size(i)
        u_resolved(start_index) = squeeze(mean(u_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        v_resolved(start_index) = squeeze(mean(v_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        w_resolved(start_index) = squeeze(mean(w_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        theta_resolved(start_index) = squeeze(mean(theta_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        start_index = start_index+1;
        end_index = end_index+1;
    end
end

u_resolved_interp = interp1(linspace(0,24,tmp),u_resolved,linspace(0,24,1728000),'spline');
v_resolved_interp = interp1(linspace(0,24,tmp),v_resolved,linspace(0,24,1728000),'spline');
w_resolved_interp = interp1(linspace(0,24,tmp),w_resolved,linspace(0,24,1728000),'spline');
theta_resolved_interp = interp1(linspace(0,24,tmp),theta_resolved,linspace(0,24,1728000),'spline');

u_sgs = u_prime(:,tower_height)' - u_resolved_interp;
v_sgs = v_prime(:,tower_height)' - v_resolved_interp;
w_sgs = w_prime(:,tower_height)' - w_resolved_interp;
theta_sgs = theta_prime(:,tower_height)' - theta_resolved_interp; 

%%
chunk = 20*60*10;
%compute terms on 30 min chunk
u_prime_sq_sgs = squeeze(mean(reshape((u_sgs.^2)',[chunk,144]),1,'omitnan'));
v_prime_sq_sgs = squeeze(mean(reshape((v_sgs.^2)',[chunk,144]),1,'omitnan'));
w_prime_sq_sgs = squeeze(mean(reshape((w_sgs.^2)',[chunk,144]),1,'omitnan'));

uw_sgs = squeeze(mean(reshape((u_sgs.*w_sgs)',[chunk,144]),1,'omitnan'));
vw_sgs = squeeze(mean(reshape((v_sgs.*w_sgs)',[chunk,144]),1,'omitnan'));
wTheta_sgs = squeeze(mean(reshape((w_sgs.*theta_sgs)',[chunk,144]),1,'omitnan'));


%Compute tke and u_star
tke_sgs = 0.5.*(u_prime_sq_sgs+v_prime_sq_sgs+...
    w_prime_sq_sgs);

u_star_sgs = (uw_sgs.^2+vw_sgs.^2).^(1/4);
%%
Pr = (uw_sgs)./(wTheta_sgs); 


%% plot
figure()
plot(linspace(0,24,144),tke_sgs(:).^(.5),'k-*')
hold on
plot(linspace(0,24,144),u_star_sgs(:),'k-o')
legend('$\sqrt{e}_{sgs}$ [m$^2$s$^{-2}$]','$u^*_{sgs}$ [ms$^{-1}$]')
xlabel("time [hrs]")
axis tight
grid on
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);
% Place second set of axes on same plot
handaxes2 = axes('Position', [0.2 0.5 0.2 0.3]);

plot(linspace(0,24,1440),SWdn,  'k--')
%set(handaxes2, 'Box', 'off')
xlabel('t (hrs)')
ylabel('SWdn [Wm$^{-2}$]')
axis tight
grid on

%% Compute pseudo flux with tower data
emiss = 0.93; 
sigma = 5.67e-8; 
z = 2.02;
z0 = 0.11e-3; %roughness length [m] for Playa (Jensen et al)
L = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.L(data_start:data_end,2:end)),linspace(0,24,1440));
T_air_tower = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.derivedT(data_start:data_end,2:4:22)),linspace(0,24,1440),'spline');
Ts = (LWup./(emiss*sigma)).^(1/4);
zeta = z./L(:,5); 

ra_shao = ra_fun(z,tke_sgs,z0,'Shao',u_star_sgs,zeta); 
ra_most = ra_fun(z,tke_sgs,z0,'MOST',u_star_sgs,zeta); 






