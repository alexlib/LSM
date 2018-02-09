%Turblence Analysis
clear all; close all; clc;
% formatting
addpath('/Users/travismorrison/Documents/Code/functions_library')
ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);
% looking at Playa on May 25th
load('/Users/travismorrison/Documents/Local_Data/MATERHORN/data/tower_data/Tower_filtered/May_2013_Playa_Filtered_N_Taylors.mat');
fprintf('Filtered Data Loaded: Filtered for Northerly Flow & u/U<0.3 \n')
load('/Users/travismorrison/Documents/Local_Data/MATERHORN/data/tower_data/Tower_filtered/May_2013_LES_filtered_data.mat')
fprintf('LES data loaded')
%% Perform Turbulence analysis
LES.sgs.uw_prime = LES.sgs.u_prime.*LES.sgs.w_prime;
LES.sgs.wT_prime = LES.sgs.w_prime.*LES.sgs.T_prime;

LES.resolved.uw_prime = LES.resolved.u_prime.*LES.sgs.w_prime;
LES.resolved.wT_prime = LES.resolved.w_prime.*LES.sgs.T_prime;



%compute derivatives
dt = 1/20;

dz_2 = 2.01-1.99;
P1 = zeros(501,1);
P2 = zeros(501,1);
for i = 1:length(LES.sgs.T(:,1))
    %fit to poly
    P1(:) = lagrangepoly([5,2.02,0.061],LES.sgs.T(i,4:6),0:.01:5);
    P2(:) = lagrangepoly([5,2.02,0.061],LES.resolved.T(i,4:6),0:.01:5);
    %take derivate as a central difference around 2.02 m. There for with
    %interpalation of 501 points we have a 2dz = .02
    %
    LES.sgs.dT_dz(i) = (P1(204) - P1(202))/dz_2;
    LES.resolved.dT_dz(i) = (P2(204) - P2(202))/dz_2;
    if mod(i,10000)==0
        index = i
    end
end

%% Compute du dz at 2m
dz_2 = 2.01-1.99;
P1 = zeros(501,1);
P2 = zeros(501,1);
for i = 1:length(LES.sgs.U(:,1))
    %fit to poly
    P1(:) = lagrangepoly([5,2.02,0.061],LES.sgs.U(i,4:6),0:.01:5);
    P2(:) = lagrangepoly([5,2.02,0.061],LES.resolved.U(i,4:6),0:.01:5);
    %take derivate as a central difference around 2.02 m. There for with
    %interpalation of 501 points we have a 2dz = .02
    %
    LES.sgs.dU_dz(i) = (P1(204) - P1(202))/dz_2;
    LES.resolved.dU_dz(i) = (P2(204) - P2(202))/dz_2;
    if mod(i,10000)==0
        index = i
    end
end

%% Compute the turbulent Pr number then Chunk it
chunk = 20*60*30;
n_30min_chunks = 531;
LES.sgs.Pr = (squeeze(mean(reshape(LES.sgs.uw_prime,[chunk,n_30min_chunks]),1,'omitnan')).*...
    squeeze(mean(reshape(LES.sgs.dT_dz,[chunk,n_30min_chunks]),1,'omitnan')))./...
    (squeeze(mean(reshape(LES.sgs.wT_prime,[chunk,n_30min_chunks]),1,'omitnan')).*...
    squeeze(mean(reshape(LES.sgs.dU_dz,[chunk,n_30min_chunks]),1,'omitnan')));

LES.resolved.Pr = (squeeze(mean(reshape(LES.resolved.uw_prime,[chunk,n_30min_chunks]),1,'omitnan')).*...
    squeeze(mean(reshape(LES.resolved.dT_dz,[chunk,n_30min_chunks]),1,'omitnan')))./...
    (squeeze(mean(reshape(LES.resolved.wT_prime,[chunk,n_30min_chunks]),1,'omitnan')).*...
    squeeze(mean(reshape(LES.resolved.dU_dz,[chunk,n_30min_chunks]),1,'omitnan')));
%% Let's examine each term first 

%velocity
figure()
subplot(3,1,1)
%
%%
chunk = 20*60*15;
n_30min_chunks = 531*2;
figure()
plot(squeeze(mean(reshape(playa_filt.Hz_20.v(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'b.-')
hold on
plot(playa_filt.min_30.w(:,tower_height),'r.-')
legend('Raw GPF $\&$ LinDetrend','LinDetrend')
ylabel('$w$')
axis tight
%%
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.u(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'b.-')
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.u(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'r.-')
grid on 
axis tight
legend('tower','resolved','SGS')
ylabel('$u$')


subplot(3,1,2)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(playa_filt.Hz_20.v(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'g.-')
hold on
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.v(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'b.-')
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.v(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'r.-')
grid on 
axis tight
%legend('tower','resolved','SGS')
ylabel('$v$')


subplot(3,1,3)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(playa_filt.Hz_20.w(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'g-*')
hold on
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.w(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'b-o')
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.w(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'r-^')
grid on 
axis tight
%legend('tower','resolved','SGS')
ylabel('$w$')
xlabel('time')

%%
subplot(3,2,2)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(playa_filt.Hz_20.u_prime(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'g*')
hold on
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.u_prime,[chunk,n_30min_chunks]),1,'omitnan')),'bo')
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.u_prime,[chunk,n_30min_chunks]),1,'omitnan')),'r^')
grid on 
axis tight
%legend('tower','resolved','SGS')
ylabel('$u^\prime$')
xlabel('time')

subplot(3,2,4)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(playa_filt.Hz_20.v_prime(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'g*')
hold on
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.v_prime,[chunk,n_30min_chunks]),1,'omitnan')),'bo')
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.v_prime,[chunk,n_30min_chunks]),1,'omitnan')),'r^')
grid on 
axis tight
%legend('tower','resolved','SGS')
ylabel('$v^\prime$')
xlabel('time')

subplot(3,2,6)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(playa_filt.Hz_20.w_prime(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'g*')
hold on
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.w_prime,[chunk,n_30min_chunks]),1,'omitnan')),'bo')
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.w_prime,[chunk,n_30min_chunks]),1,'omitnan')),'r^')
grid on 
axis tight
%legend('tower','resolved','SGS')
ylabel('$w^\prime$')
xlabel('time')

%% look at each term of the Pr sgs to see why it goes negative
tower_height = 5;
figure()
subplot(4,1,1)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape((LES.resolved.uw_prime),[chunk,n_30min_chunks]),1,'omitnan')),'k*');
ylabel('$u^\prime w^\prime$')
grid on
subplot(4,1,2)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.dT_dz,[chunk,n_30min_chunks]),1,'omitnan')),'k^');
ylabel('$\partial T / \partial z$')
grid on
subplot(4,1,3)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.wT_prime,[chunk,n_30min_chunks]),1,'omitnan')),'ko');
ylabel('$w^\prime \theta^\prime$')
grid on
subplot(4,1,4)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.dU_dz,[chunk,n_30min_chunks]),1,'omitnan')),'k>');
ylabel('$\partial U / \partial z$')
grid on
%plot(linspace(0,24,144),squeeze(mean(reshape(w_prime_sgs,[chunk,144]),1,'omitnan')));
% legend('$u^\prime w^\prime$','$\partial \theta / \partial z$',...
%     '$w^\prime \theta^\prime$','$\partial u / \partial z$','$w^\prime$')
grid on
%plot(linspace(0,24,2),[0 0],'k-')
axis tight
%%

figure()
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape((-LES.sgs.uw_prime./LES.sgs.dU_dz),[chunk,n_30min_chunks]),1,'omitnan')));
hold on
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(-LES.sgs.wT_prime./LES.sgs.dT_dz,[chunk,n_30min_chunks]),1,'omitnan')));
plot(datetime(datevec(playa_filt.min_30.t)),LES.sgs.Pr,'ko')
legend('$k_m = \frac{-u^\prime w^\prime}{\partial u / \partial z}$',...
    '$k_h = \frac{-w^\prime T^\prime}{\partial \theta / \partial z}$',...
    '$Pr = \frac{k_m}{k_h}$')
grid on
%plot(linspace(0,24,2),[0 0],'k-')
axis tight



%%

g = 9.8; 
%Compute Rf #
LES.sgs.Rf = squeeze(mean(reshape((g.*LES.sgs.wT_prime)./(LES.sgs.T(:,tower_height)'...
    .*LES.sgs.uw_prime.*LES.sgs.dU_dz),[chunk,n_30min_chunks]),1,'omitnan'));

%% Compute tke and u_star


%compute terms on 30 min chunk
u_prime_sq_sgs = squeeze(mean(reshape((LES.sgs.u_prime.^2)',[chunk,n_30min_chunks]),1,'omitnan'));
v_prime_sq_sgs = squeeze(mean(reshape((LES.sgs.v_prime.^2)',[chunk,n_30min_chunks]),1,'omitnan'));
w_prime_sq_sgs = squeeze(mean(reshape((LES.sgs.w_prime.^2)',[chunk,n_30min_chunks]),1,'omitnan'));

LES.sgs.tke = 0.5.*(u_prime_sq_sgs+v_prime_sq_sgs+...
    w_prime_sq_sgs);
%ustar vars
uw_sgs = squeeze(mean(reshape((LES.sgs.u_prime.*LES.sgs.w_prime)',[chunk,n_30min_chunks]),1,'omitnan'));
vw_sgs = squeeze(mean(reshape((LES.sgs.v_prime.*LES.sgs.w_prime)',[chunk,n_30min_chunks]),1,'omitnan'));

LES.sgs.u_star = (uw_sgs.^2+vw_sgs.^2).^(1/4);


%%

%% Look at heat flux of resolved vs sgs vs tower
figure()
plot(squeeze(mean(reshape(playa_filt.Hz_20.w_prime(:,5).*playa_filt.Hz_20.T_prime(:,5),[chunk,n_30min_chunks]),1,'omitnan')));
hold on 
plot(squeeze(mean(reshape(LES.resolved.wT_prime,[chunk,n_30min_chunks]),1,'omitnan')));
plot(squeeze(mean(reshape(LES.sgs.wT_prime,[chunk,n_30min_chunks]),1,'omitnan')));
plot(squeeze(mean(reshape(LES.resolved.wT_prime+LES.sgs.wT_prime,[chunk,n_30min_chunks]),1,'omitnan')));
grid on 
axis tight
legend('tower','resolved','SGs','resolved+sgs')
ylabel('$w^\prime T^\prime$')
xlabel('time ')
%% plot coeffs vs time
figure()
%plot(datetime(datevec(playa_filt.min_30.t)),((Ck.*LES.sgs.tke.^(.5))./LES.sgs.Pr),'ko')
hold on 
%plot(datetime(datevec(playa_filt.min_30.t)),(kappa.*LES.sgs.u_star),'k*')
plot(datetime(datevec(playa_filt.min_30.t)),LES.resolved.Pr,'k*')
grid on 
axis tight
%legend('$C_k \sqrt{e_{sgs}} Pr^{-1}$','$\kappa u^*_{sgs}$')
ylabel('$Pr_{resolved}$')
xlabel('time ')
%% plot Pr vs Rf
figure()
s(LES.sgs.Rf, LES.sgs.Pr,'ko')
% hold on 
% plot([-10e5 10e5], [0 0],'k-')
% plot([0 0],[-10e5 10e5], 'k-')
% xlim([-10 10])
% ylim([-10 10])
ylabel('$Pr_{sgs}$')
xlabel('$R_{f}$')
grid on
axis tight
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);



%% Plot kappa ustar against tke/Pr
kappa =0.4;
Ck= 0.15;
min = 0; 
max = 1.5;
figure()
plot(((Ck.*LES.sgs.tke.^(.5))./.3),(kappa.*LES.sgs.u_star),'ko')
hold on
plot(linspace(-10000,10000,2),linspace(-10000,10000,2),'k--')
title('Pr  = 0.3','interpreter','latex')
ylim([min  max])
xlim([min max])
xlabel("")
%axis tight
grid on

ylabel('$\kappa u^*_{sgs}$')
xlabel('$C_k \sqrt{e_{sgs}} Pr^{-1}$')