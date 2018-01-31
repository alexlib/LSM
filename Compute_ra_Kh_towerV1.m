%% Compute r_h and Kh for the tower as a function of height
%following Shao 2013
clear all; close all; clc;
% formatting

%addpath('/Users/travismorrison/Documents/Code/functions_library')
ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);
%%
%load tower data: 30 min data
data_start = 1072;%25 May 0000 UTC
data_end = data_start + (2*24);%26 May 0000 UTC
load('.\Materhorn_data\playaSpring30minLinDetUTESpac3.mat')
U_interp = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.spdAndDir(data_start:data_end,3:3:18)),linspace(0,24,1440),'spline');
U = rearrangeHeights(playaSpring.spdAndDir(data_start:data_end,3:3:18));

tke = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.tke(data_start:data_end,2:end)),linspace(0,24,1440));
L = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.L(data_start:data_end,2:end)),linspace(0,24,1440));
T_air_tower = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.derivedT(data_start:data_end,2:4:22)),linspace(0,24,1440),'spline');
%%


%%
%Examine velocity based on different averaging
figure()
plot(linspace(0,24.5,1440),U_interp(:,1))
hold on
plot(linspace(0,24.5,49),U(:,1))
plot(linspace(0,24.5,24),U_1hr(:,1),'k*')
legend('UtesPac','LSM','1hr 3pt interp')
xlabel('time (hrs)')
ylabel('U$_{25m}$ (ms$^{-1}$)')
axis tight

%look at tke values 
figure()
plot(linspace(0,24.5,1440),tke(:,1))

legend('30min 3pt interp')
xlabel('time (hrs)')
ylabel('TKE$_{25m}$ (ms$^{-1}$)')
axis tight
%% apply talyor's forzen turbo hypo
time_chunk = 60*60; %[sec] currently 1hr chunks
x = (tke.*time_chunk);
x_1hr = reshape(x,[(size(tke,1)/24),24,6]);
x_1hr_mean = squeeze(mean(x_1hr,1));

% look at distance of chunks 
figure()
plot(linspace(0,24,24),x_1hr_mean(:,1),'k*')
axis tight
legend('30min 3pt interp')
xlabel('time (hrs)')
ylabel('x$_{25m}$ (m)')
axis tight

%Apply filter on LES scale (sliding average)
LES_res = 10; %[m]
tke_1hr = reshape(tke,[(size(tke,1)/24),24,6]);
nx = floor(x_1hr_mean/LES_res); %Think LES domain nx, number of x grid points (coarse grid)
% 
% for z = 1:6
%     for 
%     end
% end
% 
% 


