%following Shao 2013
clear all; close all; clc;
% formatting
addpath('/Users/travismorrison/Documents/Code/functions_library')
ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);
%% looking at Playa on May 25th
load('/Users/travismorrison/Documents/Local_Data/MATERHORN/data/tower_data/Playa_tower_raw/PlayaSpring_raw_2013_05_24.mat')
start_index = 432000; %May 24 0600 UTC or 0000MST
end_index = start_index +(20*60*60*24); %May 25 0800 UTC or 0200MST
% load radiation data: 5 min data
data_start = 6696; %25 May 0600 UTC
data_end = data_start + (12*24); %26 May 0800 UTC
load('./Materhorn_data/MATERHORN_Rad_data.mat');
SWdn = interp1(linspace(0,24,289),rad_data(data_start:data_end,6),linspace(0,24,1440),'spline');
%%
u_prime = rearrangeHeights(rawFlux.uPrime(start_index:end_index-1,:));
v_prime = rearrangeHeights(rawFlux.vPrime(start_index:end_index-1,:));
w_prime = rearrangeHeights(rawFlux.wPrime(start_index:end_index-1,:));
chunk = 20*60*10;

%compute terms on 30 min chunk
u_prime_sq = squeeze(mean(reshape((u_prime.^2)',[chunk,144,6]),1,'omitnan'));
v_prime_sq = squeeze(mean(reshape((v_prime.^2)',[chunk,144,6]),1,'omitnan'));
w_prime_sq = squeeze(mean(reshape((w_prime.^2)',[chunk,144,6]),1,'omitnan'));

uw_prime = squeeze(mean(reshape((u_prime.*w_prime)',[chunk,144,6]),1,'omitnan'));
vw_prime = squeeze(mean(reshape((v_prime.*w_prime)',[chunk,144,6]),1,'omitnan'));


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

%% take 30 min chunk average

