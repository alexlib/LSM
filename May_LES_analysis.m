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

%%


wT_sgs = squeeze(mean(reshape(playa_filt.Hz_20.wT_prime(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan'))-...
    (squeeze(mean(reshape(LES.resolved.w(:,tower_height).*LES.resolved.T(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan'))-...
    squeeze(mean(reshape(LES.resolved.w(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')).* ...
    squeeze(mean(reshape(LES.resolved.T(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')));




%% compute derivatives
dt = 1/20;

dz_2 = 2.01-1.99;
P1 = zeros(501,1);
P2 = zeros(501,1);
leng = length(LES.sgs.T);
LES.sgs.dT_dz = zeros(leng,1);
LES.resolved.dT_dz = zeros(leng,1);
for i = 1:length(LES.sgs.T(:,1))
    %fit to poly
    P1(:) = lagrangepoly([5,2.02,0.61],LES.sgs.T(i,4:6),0:.01:5);
    P2(:) = lagrangepoly([5,2.02,0.61],LES.resolved.T(i,4:6),0:.01:5);
    %take derivate as a central difference around 2.02 m. There for with
    %interpalation of 501 points we have a 2dz = .02
    %
    LES.sgs.dT_dz(i) = (P1(204) - P1(202))/dz_2;
    LES.resolved.dT_dz(i) = ( P2(204) -  P2(202))/dz_2;
    if mod(i,10000)==0
        index = i
    end
end
fprintf('dTdz complete \n')

% Compute du dz at 2m
dz_2 = 2.01-1.99;
P1 = zeros(501,1);
P2 = zeros(501,1);
LES.sgs.dU_dz = zeros(leng,1);
LES.resolved.dU_dz = zeros(leng,1);
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
fprintf('dUdz complete \n')
%% Compute the turbulent Pr number then Chunk it
chunk = 20*60*30;
n_30min_chunks = 580;
LES.sgs.Pr = (squeeze(mean(reshape(LES.sgs.uw_prime,[chunk,n_30min_chunks]),1,'omitnan')).*...
    squeeze(mean(reshape(LES.sgs.dT_dz,[chunk,n_30min_chunks]),1,'omitnan')))./...
    (squeeze(mean(reshape(LES.sgs.wT_prime,[chunk,n_30min_chunks]),1,'omitnan')).*...
    squeeze(mean(reshape(LES.sgs.dU_dz,[chunk,n_30min_chunks]),1,'omitnan')));

LES.resolved.Pr = (squeeze(mean(reshape(LES.resolved.uw_prime,[chunk,n_30min_chunks]),1,'omitnan')).*...
    squeeze(mean(reshape(LES.resolved.dT_dz,[chunk,n_30min_chunks]),1,'omitnan')))./...
    (squeeze(mean(reshape(LES.resolved.wT_prime,[chunk,n_30min_chunks]),1,'omitnan')).*...
    squeeze(mean(reshape(LES.resolved.dU_dz,[chunk,n_30min_chunks]),1,'omitnan')));
%% Let's examine each term first 
% double check u and T
chunk = 20*60*30;
n_30min_chunks = 580;
tower_height = 5;

figure()
subplot(2,1,1)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.u(:,tower_height),...
    [chunk,n_30min_chunks]),1,'omitnan')),'b.-')
hold on
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.u(:,tower_height),...
    [chunk,n_30min_chunks]),1,'omitnan')),'r.-')
plot(datetime(datevec(playa_filt.min_30.t)),playa_filt.min_30.u(:,tower_height),'k-');
ylabel('u [m/s]')
legend('SGS','Resolved','Total')
axis tight

subplot(2,1,2)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.T(:,tower_height),...
    [chunk,n_30min_chunks]),1,'omitnan')),'b.-')
hold on
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.T(:,tower_height),...
    [chunk,n_30min_chunks]),1,'omitnan')),'r.-')
plot(datetime(datevec(playa_filt.min_30.t)),playa_filt.min_30.T(:,tower_height),'k-');
ylabel('T [K]')
axis tight
%% Examine the u'w' and w'T'
figure()
subplot(2,1,1)
l1 = plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.uw_prime(:),...
    [chunk,n_30min_chunks]),1,'omitnan')),'b*');
hold on
l2 = plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.uw_prime(:),...
    [chunk,n_30min_chunks]),1,'omitnan')),'r*');
l3 = plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(playa_filt.Hz_20.u_prime(:,tower_height).*playa_filt.Hz_20.w_prime(:,tower_height), [chunk,n_30min_chunks]),1,'omitnan')),'k*');
ylabel('$\overline{u^\prime w^\prime}$ $(m^{2}s^{-2})$');
legend('SGS','Resolved','Total')
axis tight
grid on
l1.MarkerSize = 4;
l2.MarkerSize = 4;
l3.MarkerSize = 4;

subplot(2,1,2)
l1 = plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.wT_prime,...
    [chunk,n_30min_chunks]),1,'omitnan')),'b*');
hold on
l2 = plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.wT_prime,...
    [chunk,n_30min_chunks]),1,'omitnan')),'r*');
l3 = plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(playa_filt.Hz_20.wT_prime(:,tower_height),[chunk,n_30min_chunks]),1,'omitnan')),'k*');
plot(datetime(datevec(playa_filt.min_30.t)),wT_sgs)

ylabel('$\overline{w^\prime T^\prime}$ $(ms^{-1}K)$')
axis tight
grid on
l1.MarkerSize = 4;
l2.MarkerSize = 4;
l3.MarkerSize = 4;


%%
figure()
subplot(4,2,1)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.uw_prime,[chunk,n_30min_chunks]),1,'omitnan')),'ko');
ylabel('$w^\prime \theta^\prime$')
grid on
title('SGS')
subplot(4,2,3)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.dT_dz,[chunk,n_30min_chunks]),1,'omitnan')),'k^');
ylabel('$\partial T / \partial z$')
grid on
subplot(4,2,5)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.wT_prime,[chunk,n_30min_chunks]),1,'omitnan')),'ko');
ylabel('$w^\prime \theta^\prime$')
grid on
subplot(4,2,7)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.sgs.dU_dz,[chunk,n_30min_chunks]),1,'omitnan')),'k>');
ylabel('$\partial U / \partial z$')
grid on
axis tight

subplot(4,2,2)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.uw_prime,[chunk,n_30min_chunks]),1,'omitnan')),'ko');
ylabel('$w^\prime \theta^\prime$')
grid on
title('Resolved')
subplot(4,2,4)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.dT_dz,[chunk,n_30min_chunks]),1,'omitnan')),'k^');
ylabel('$\partial T / \partial z$')
grid on
subplot(4,2,6)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.wT_prime,[chunk,n_30min_chunks]),1,'omitnan')),'ko');
ylabel('$w^\prime \theta^\prime$')
grid on
subplot(4,2,8)
plot(datetime(datevec(playa_filt.min_30.t)),squeeze(mean(reshape(LES.resolved.dU_dz,[chunk,n_30min_chunks]),1,'omitnan')),'k>');
ylabel('$\partial U / \partial z$')
grid on
axis tight


%% look at T-profile with time - Make Movie
Temp1 = squeeze(mean(reshape(LES.sgs.T, [chunk,n_30min_chunks,6]),1,'omitnan'));
Temp2 = squeeze(mean(reshape(LES.resolved.T, [chunk,n_30min_chunks,6]),1,'omitnan'));
Temp3 = squeeze(mean(reshape(playa_filt.Hz_20.T, [chunk,n_30min_chunks,6]),1,'omitnan'));
for i = 1:580
    P(i,:)  = lagrangepoly([5,2.02,0.61],Temp3(i,4:6),0:.01:5);
    P2(i,:)  = lagrangepoly([5,2.02,0.61],Temp1(i,4:6),0:.01:5);
    P3(i,:)  = lagrangepoly([5,2.02,0.61],Temp2(i,4:6),0:.01:5);
end

% Make movie 
% Prepare the new file.
%vidObj = VideoWriter('60min_decomp_filtered@1Hz_sampled@0_25seconds_video.avi');
%open(vidObj); 


%fig = figure('rend','painters','pos',[10 10 600 600]);
fig = figure('pos',[10 10 600 600]);
F(:)=getframe(fig);
close(fig);
fig = figure('pos',[10 10 600 600]);
indexmovie=1:5:12000;
for i = 271:n_30min_chunks
    subplot(3,1,1)
    l1 = plot(squeeze(P(i,:)),linspace(0,5,501),'-');
    l1.Color=[.7 .7 .7];
    hold on 
    plot(Temp3(i,4:6),[5,2.02,0.61],'k*')
    ylabel('$z_{raw}$')
     hold off
     
    subplot(3,1,2)
    l1 = plot(squeeze(P2(i,:)),linspace(0,5,501),'-');
    l1.Color=[.7 .7 .7];
    hold on
    plot(Temp1(i,4:6),[5,2.02,0.61],'k*')
    ylabel('$z_{sgs}$')
     hold off
     
    subplot(3,1,3)
    l1 = plot(squeeze(P3(i,:)),linspace(0,5,501),'-');
    l1.Color=[.7 .7 .7];
    hold on
    plot(Temp2(i,4:6),[5,2.02,0.61],'k*')
    ylabel('$z_{resolved}$')
    xlabel('T [k]')
    pause(0.5)
     hold off
     
     F(i) = getframe(fig);
      
end

writerObj=VideoWriter('T_profile_total_sgs_resolved.avi','Uncompressed AVI');
writerObj.FrameRate = 4;
open(writerObj)
writeVideo(writerObj,F)
close(writerObj)

%%

g = 9.8; 
%Compute Rf #
LES.sgs.Rf = squeeze(mean(reshape((g.*LES.sgs.wT_prime)./(LES.sgs.T(:,tower_height)'...
    .*LES.sgs.uw_prime.*LES.sgs.dU_dz'),[chunk,n_30min_chunks]),1,'omitnan'));
%Compute Rf #
LES.resolved.Rf = squeeze(mean(reshape((g.*LES.resolved.wT_prime)./(LES.resolved.T(:,tower_height)'...
    .*LES.resolved.uw_prime.*LES.resolved.dU_dz'),[chunk,n_30min_chunks]),1,'omitnan'));
%% Plot Pr & Rf with time
figure()
subplot(2,1,1)
l1 = semilogy(datetime(datevec(playa_filt.min_30.t)),LES.sgs.Pr,'k*');
hold on
l2 = semilogy(datetime(datevec(playa_filt.min_30.t)),LES.resolved.Pr,'b*');
axis tight
grid on
legend('SGS','Resolved')
ylabel('Pr')

subplot(2,1,2)
l1 = semilogy(datetime(datevec(playa_filt.min_30.t)),LES.sgs.Rf,'k*');
hold on
l2 = semilogy(datetime(datevec(playa_filt.min_30.t)),LES.resolved.Rf,'b*');
axis tight
grid on
ylabel('$R_f$')

%% Pr vs Rf and Pr vs L
figure()
subplot(2,2,1)

loglog(LES.sgs.Rf,LES.sgs.Pr,'k*')
grid on
axis tight
xlabel('$R_f$')
ylabel('$Pr$')
title('SGS')
subplot(2,2,3)
loglog(playa_filt.min_30.L,LES.sgs.Pr,'k*')
grid on
axis tight
xlabel('$L$')
ylabel('$Pr$')

subplot(2,2,2)
loglog(LES.resolved.Rf,LES.resolved.Pr,'k*')
grid on
axis tight
xlabel('$R_f$')
%ylabel('$Pr$')
title('Resolved')
subplot(2,2,4)
loglog(playa_filt.min_30.L,LES.resolved.Pr,'k*')
grid on
axis tight
xlabel('$L$')
%ylabel('$Pr$')



%% Compute tke and u_star - These terms need to be calcualted before... 


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