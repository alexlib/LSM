%Prepare_data_Test
%load data
clear; clc;
% settings header
addpath('/Users/travismorrison/Documents/Code/functions_library')
ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);
%Import data which has been properly orientated by UTESPAC (LPF on 30 minute chunks and GPF), filtered for
%Northerly flow and meets taylors hypothesis requirment
load('/Users/travismorrison/Documents/Local_Data/MATERHORN/data/tower_data/Tower_filtered/May_2013_Playa_Filtered_N_Taylors.mat');
fprintf('Filtered Data Loaded: Filtered for Northerly Flow & u/U<0.3 \n')

% reshape data into time series from 30min chunks
z = [25.5,19.4,10,5,2.02,0.61];
tower_height = 5; %2-m height index
playa_filt.Hz_20.w = reshape(playa_filt.Hz_20.w, [36000*531,6]);
playa_filt.Hz_20.u_prime = reshape(playa_filt.Hz_20.u_prime, [36000*531,6]);
playa_filt.Hz_20.v_prime = reshape(playa_filt.Hz_20.v_prime, [36000*531,6]);
playa_filt.Hz_20.T = reshape(playa_filt.Hz_20.T, [36000*531,6]);
for n_chunks = 1:531
for i = 1:length(z)
    playa_filt.Hz_20.Th_prime(:,n_chunks,i) = detrend(squeeze(playa_filt.Hz_20.Th(:,n_chunks,i)));
end
end
playa_filt.Hz_20.Th = reshape(playa_filt.Hz_20.Th,[36000*531,6]);
playa_filt.Hz_20.Th_prime = reshape(playa_filt.Hz_20.Th_prime,[36000*531,6]);
playa_filt.Hz_20.u = reshape(playa_filt.Hz_20.u, [36000*531,6]);
playa_filt.Hz_20.v = reshape(playa_filt.Hz_20.v, [36000*531,6]);
playa_filt.Hz_20.w_prime = reshape(playa_filt.Hz_20.w_prime, [36000*531,6]);
playa_filt.Hz_20.T_prime = reshape(playa_filt.Hz_20.T_prime, [36000*531,6]);



% load radiation data: 5 min data for estimate of surface temperature
load('./Materhorn_data/MATERHORN_Rad_data.mat');
% From looking at data set, the entire day of May 10th was a solid day...
% Lets check the Obukhov Length for validity of MOST
%enter a day which has a continous chunk of afternoon data
%
day_of_interest = 'May14';
switch day_of_interest
    case 'May14'%May 14th
        %start_i = 91;
        start_i = 100;
        end_i = 104;
        rad_index = [3953:1:3965];
        %rad_index = [3881:1:3965];
    case 'May15'%May 15th
        start_i = 107;
        end_i = 137;
        rad_index = [4067:1:4247];
    case 'May17'%May 17th
        start_i = 180;
        end_i = 209;
    case 'May20'%Mary 20th
        start_i = 247;
        end_i = 275;
        rad_index = [5561:5729];
    case 'May31'%Mary 31th
        start_i = 474;
        end_i = 483;
end
start_i20hz = start_i*20*30*60;
end_i20hz = end_i*20*30*60-1;

%load radition data
LWup = interp1(linspace(0,24,size(rad_index,2)),rad_data(rad_index,9),linspace(0,24,(end_i20hz-start_i20hz+1)),'spline');
SWdn = interp1(linspace(0,24,size(rad_index,2)),rad_data(rad_index,6),linspace(0,24,(end_i20hz-start_i20hz+1)),'spline');
emiss = .97;
sb= 5.67E-8;
playa_filt.Hz_20.T_sur = (LWup./(emiss*sb)).^(1/4)-273.15;
playa_filt.min_30.T_sur = squeeze(mean(reshape(playa_filt.Hz_20.T_sur,[20*30*60,(end_i-start_i)]),1,'omitnan'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%look at L
figure()
subplot(3,2,1)
plot(datetime(datevec(playa_filt.min_30.t(start_i:end_i))),z(tower_height)./...
    playa_filt.min_30.L(start_i:end_i,tower_height),'k-o')
grid on
ylabel('$\zeta$')
title(day_of_interest)

%Subsdiance? 
subplot(3,2,2)
plot(linspace(0,(end_i-start_i)*30,(end_i20hz-start_i20hz+1)),playa_filt.Hz_20.w(start_i20hz:end_i20hz,tower_height)','k-')
grid on
ylabel('w')
xlabel('time (minutes)')
axis tight

%constant flux layer? 
wTh_prime = squeeze(mean(reshape(playa_filt.Hz_20.w_prime(start_i20hz:end_i20hz,:).*...
    playa_filt.Hz_20.Th_prime(start_i20hz:end_i20hz,:),[20*60*30,(end_i-start_i),6]),1,'omitnan'));
uw_prime = squeeze(mean(reshape(playa_filt.Hz_20.u_prime(start_i20hz:end_i20hz,:).*...
    playa_filt.Hz_20.w_prime(start_i20hz:end_i20hz,:),[20*60*30,(end_i-start_i),6]),1,'omitnan'));

subplot(3,2,3)
plot(wTh_prime,z,'k-o')
xlabel('$\overline{w^\prime \theta^\prime}$')
ylabel('z (m)')
grid on

subplot(3,2,4)
plot(uw_prime,z,'k-o')
xlabel('$\overline{u^\prime w^\prime}$')
ylabel('z (m)')
grid on

subplot(3,2,5)
plot(datetime(datevec(playa_filt.min_30.t(start_i:end_i-1))),playa_filt.min_30.T_sur)
grid on 
axis tight
ylabel('$T_s$')

subplot(3,2,6)
plot(SWdn)
grid on 
axis tight
ylabel('$S_{\downarrow}$')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute MRD on 1Hr chunks for wT & uw 
M = floor(log2(20*60*60));
n_chunks = floor((end_i-start_i)./2);
mrd_start_i = start_i20hz;
mrd_end_i = start_i20hz+(20*60*60);

for  h = 1:length(z)
    for n = 1:n_chunks
        uw_MRD(:,h,n) = compute_MRD(playa_filt.Hz_20.u_prime(mrd_start_i:mrd_end_i,h),...
            playa_filt.Hz_20.w_prime(mrd_start_i:mrd_end_i,h),M,0);
        wT_MRD(:,h,n) = compute_MRD(playa_filt.Hz_20.w_prime(mrd_start_i:mrd_end_i,h),...
            playa_filt.Hz_20.T_prime(mrd_start_i:mrd_end_i,h),M,0);
        mrd_start_i = mrd_end_i+1;
        mrd_end_i = mrd_end_i+(20*60*60);
    end
    mrd_start_i = start_i20hz;
    mrd_end_i = start_i20hz+(20*60*60);
end

%take the mean over all the 1 hr chunks
uw_MRD_mean = squeeze(mean(uw_MRD,3,'omitnan'));
wT_MRD_mean = squeeze(mean(wT_MRD,3,'omitnan'));


%take the standard deviation over all the periods
uw_MRD_std = std(uw_MRD,0,3,'omitnan');
wT_MRD_std = std(wT_MRD,0,3,'omitnan');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the MRD
x_axis = (2.^[1:M]/(20));
xmin = x_axis(1);
xmax = x_axis(end);

figure()
subplot(1,2,1)
errorbar(x_axis,uw_MRD_mean(:,1),uw_MRD_std(:,1),'k-o')
hold on
errorbar(x_axis,uw_MRD_mean(:,2),uw_MRD_std(:,2),'b-o')
errorbar(x_axis,uw_MRD_mean(:,3),uw_MRD_std(:,3),'y-o')
errorbar(x_axis,uw_MRD_mean(:,4),uw_MRD_std(:,4),'c-o')
errorbar(x_axis,uw_MRD_mean(:,5),uw_MRD_std(:,5),'r-o')
errorbar(x_axis,uw_MRD_mean(:,6),uw_MRD_std(:,6),'g-o')
set(gca,'xscale','log')
%set(gca,'yscale','log')
grid on
axis tight
title('$u^\prime w^\prime$')
legend('25-m','19','10','5','2.02','0.61','Location','southwest')


subplot(1,2,2)
errorbar(x_axis,wT_MRD_mean(:,1),wT_MRD_std(:,1),'k-o')
hold on
errorbar(x_axis,wT_MRD_mean(:,2),wT_MRD_std(:,2),'b-o')
errorbar(x_axis,wT_MRD_mean(:,3),wT_MRD_std(:,3),'y-o')
errorbar(x_axis,wT_MRD_mean(:,4),wT_MRD_std(:,4),'c-o')
errorbar(x_axis,wT_MRD_mean(:,5),wT_MRD_std(:,5),'r-o')
errorbar(x_axis,wT_MRD_mean(:,6),wT_MRD_std(:,6),'g-o')
set(gca,'xscale','log')
%set(gca,'yscale','log')
grid on 
axis tight
title('$w^\prime T^\prime$')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute MOST and Shao flux on 30 minute averages
Pr = 0.3; %Turbulent Pr number from Shao et al 
z0 = 0.11e-3; %roughness length [m] for Playa (Jensen et al)
for h = 1:length(z)
playa_filt.Hz_20.u_star(:,h) = ((playa_filt.Hz_20.u_prime(start_i20hz:end_i20hz,h).*...
    playa_filt.Hz_20.w_prime(start_i20hz:end_i20hz,h)).^2).^(1/4);
end
playa_filt.min_30.u_star = squeeze(mean(reshape(playa_filt.Hz_20.u_star,...
    [20*30*60,(end_i-start_i),6]),1,'omitnan'));
n_chunks = floor((end_i-start_i));
flux_start_i = start_i20hz;
flux_end_i = start_i20hz+(20*60*30);

for  h = 1:length(z)
    
    MOST_ra(:,h) = ra_fun(z(h),playa_filt.min_30.tke(start_i:end_i-1,h),Pr,z0,'MOST',...
        playa_filt.min_30.u_star(:,h),z(h)./playa_filt.min_30.L(start_i:end_i,h));
    Shao_ra(:,h) = ra_fun(z(h),playa_filt.min_30.tke(start_i:end_i-1,h),Pr,z0,'Shao',...
        playa_filt.min_30.u_star(:,h),z(h)./playa_filt.min_30.L(start_i:end_i,h));
 
end
%
rho_Cp_air = 1.1614.*1.007e3;
MOST_H = -rho_Cp_air.*((playa_filt.min_30.T(start_i:end_i-1,tower_height) - playa_filt.min_30.T_sur')./MOST_ra(:,tower_height));
Shao_H = -rho_Cp_air.*((playa_filt.min_30.T(start_i:end_i-1,tower_height) - playa_filt.min_30.T_sur')./Shao_ra(:,tower_height));
playa_filt.min_30.H = rho_Cp_air.*flip(wTh_prime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure (3) Fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(datetime(datevec(playa_filt.min_30.t(start_i:end_i-1))),MOST_H,'r-o')
hold on
plot(datetime(datevec(playa_filt.min_30.t(start_i:end_i-1))),Shao_H,'b-o')
plot(datetime(datevec(playa_filt.min_30.t(start_i:end_i-1))),playa_filt.min_30.H(:,tower_height),'k-o')
legend('MOST','Shao','Tower')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

