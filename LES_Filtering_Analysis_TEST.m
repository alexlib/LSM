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
%load data
load('/Users/travismorrison/Documents/Local_Data/MATERHORN/data/tower_data/Tower_filtered/May_2013_Playa_Filtered_N_Taylors.mat');
fprintf('Filtered Data Loaded: Filtered for Northerly Flow & u/U<0.3 \n')

% reshape data into time series
tower_height = 5;
playa_filt.Hz_20.w = reshape(playa_filt.Hz_20.w, [36000*531,6]);
playa_filt.Hz_20.u_prime = reshape(playa_filt.Hz_20.u_prime, [36000*531,6]);
playa_filt.Hz_20.v_prime = reshape(playa_filt.Hz_20.v_prime, [36000*531,6]);
playa_filt.Hz_20.T = reshape(playa_filt.Hz_20.T, [36000*531,6]);
playa_filt.Hz_20.Th = reshape(playa_filt.Hz_20.Th,[36000*531,6]);
playa_filt.Hz_20.u = reshape(playa_filt.Hz_20.u, [36000*531,6]);
playa_filt.Hz_20.v = reshape(playa_filt.Hz_20.v, [36000*531,6]);
playa_filt.Hz_20.w_prime = reshape(playa_filt.Hz_20.w_prime, [36000*531,6]);
playa_filt.Hz_20.T_prime = reshape(playa_filt.Hz_20.T_prime, [36000*531,6]);


% load radiation data: 5 min data for estimate of surface temperature
data_index = [3155:1:3160 3173:3184]; 
load('./Materhorn_data/MATERHORN_Rad_data.mat');
LWup = interp1(linspace(0,24,18),rad_data(data_index,9),linspace(0,24,108000),'spline');
emiss = .97;
sb= 5.67E-8;
playa_filt.Hz_20.T_sur = (LWup./(emiss*sb)).^(1/4)-273.15;
playa_filt.min_30.T_sur = squeeze(mean(reshape(playa_filt.Hz_20.T_sur,[36000,3]),1,'omitnan'));
clear rad_data

%Apply Taylors frozen Turbo Hypothesis
chunk_30min_start = 77;
chunks = 3;
chunks_30min_end = chunk_30min_start+chunks;
npoints = 108000;
n_start = 2664000;
n_end = n_start+npoints-1;

freq = 20;
Resolved.T = zeros(npoints,6);
for k = 1:6
   [tmp1, filter_size]= Taylors_filter( playa_filt.min_30.U(chunk_30min_start:chunks_30min_end,k),...
       playa_filt.Hz_20.T(n_start:n_end,k), freq );
   [tmp2, filter_size]= Taylors_filter( playa_filt.min_30.U(chunk_30min_start:chunks_30min_end,k),...
       playa_filt.Hz_20.u(n_start:n_end,k), freq );
    Resolved.T(:,k) = interp1(linspace(0,1,(npoints - (sum(floor(filter_size(1:chunks))))-chunks)),tmp1',linspace(0,1,npoints));
    Resolved.u(:,k) = interp1(linspace(0,1,(npoints - (sum(floor(filter_size(1:chunks))))-chunks)),tmp2',linspace(0,1,npoints));
    SFS.T(:,k) = playa_filt.Hz_20.T(n_start:n_end,k) - Resolved.T(:,k);
    clear tmp1 tmp2;
end
[ Resolved.w ] = Taylors_filter( playa_filt.min_30.U(chunk_30min_start:chunks_30min_end,tower_height),...
   playa_filt.Hz_20.w(n_start:n_end,tower_height), freq );

[ Resolved.u ] = Taylors_filter( playa_filt.min_30.U(chunk_30min_start:chunks_30min_end,tower_height),...
   playa_filt.Hz_20.u(n_start:n_end,tower_height), freq );

[ Resolved.v ] = Taylors_filter( playa_filt.min_30.U(chunk_30min_start:chunks_30min_end,tower_height),...
   playa_filt.Hz_20.v(n_start:n_end,tower_height), freq );

[ Resolved.Th ] = Taylors_filter( playa_filt.min_30.U(chunk_30min_start:chunks_30min_end,tower_height),...
   playa_filt.Hz_20.Th(n_start:n_end,tower_height), freq );

[ Resolved.wTh ] = Taylors_filter( playa_filt.min_30.U(chunk_30min_start:chunks_30min_end,tower_height),...
   playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.Th(n_start:n_end,tower_height), freq );

[ Resolved.wT ] = Taylors_filter( playa_filt.min_30.U(chunk_30min_start:chunks_30min_end,tower_height),...
   playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end,tower_height), freq );

[ Resolved.wT_prime, filter_size,psedo_space  ] = Taylors_filter( playa_filt.min_30.U(chunk_30min_start:chunks_30min_end,tower_height),...
   playa_filt.Hz_20.w_prime(n_start:n_end,tower_height).*playa_filt.Hz_20.T_prime(n_start:n_end,tower_height), freq );

%Reinterpolate back to orginal size (doesn't really effect data, just makes it easier to work with)
Resolved.w =interp1(linspace(0,1,(npoints - (sum(floor(filter_size(1:chunks))))-chunks)),Resolved.w,linspace(0,1,npoints));
Resolved.u = interp1(linspace(0,1,(npoints - (sum(floor(filter_size(1:chunks))))-chunks)),Resolved.u,linspace(0,1,npoints));
Resolved.v = interp1(linspace(0,1,(npoints - (sum(floor(filter_size(1:chunks))))-chunks)),Resolved.v,linspace(0,1,npoints));
Resolved.Th = interp1(linspace(0,1,(npoints - (sum(floor(filter_size(1:chunks))))-chunks)),Resolved.Th,linspace(0,1,npoints));
Resolved.wTh = interp1(linspace(0,1,(npoints - (sum(floor(filter_size(1:chunks))))-chunks)),Resolved.wTh,linspace(0,1,npoints));
Resolved.wT = interp1(linspace(0,1,(npoints - (sum(floor(filter_size(1:chunks))))-chunks)),Resolved.wT,linspace(0,1,npoints));

%compute SFS Scale
SFS.w = playa_filt.Hz_20.w(n_start:n_end,tower_height)' - Resolved.w;
SFS.u = playa_filt.Hz_20.u(n_start:n_end,tower_height)' - Resolved.u;
SFS.v = playa_filt.Hz_20.v(n_start:n_end,tower_height)' - Resolved.v;
SFS.Th = playa_filt.Hz_20.Th(n_start:n_end,tower_height)' - Resolved.Th;
SFS.wTh = playa_filt.Hz_20.w(n_start:n_end,tower_height)'.*playa_filt.Hz_20.Th(n_start:n_end,tower_height)' - Resolved.wTh;
SFS.wT = playa_filt.Hz_20.w(n_start:n_end,tower_height)'.*playa_filt.Hz_20.T(n_start:n_end,tower_height)' - Resolved.wT;

%% Examine the velocity and T filtered variables (Check,/makes sense)
% from here it is pertant to not that the mean of the resolved field and
% the orginal data is not 0 for all variables, like u, v and w. However the
% SFS field displays appropaite behavior (~ around 0).
figure()
subplot(4,1,1)
plot(psedo_space,SFS.u,'r-')
hold on 
plot(psedo_space,playa_filt.Hz_20.u(n_start:n_end,tower_height),'k-')
plot(psedo_space,Resolved.u,'b-')
legend('SFS','Total','Resolved')
grid on
ylabel('u [ms$^{-1}$]')
title('Filtered Velocities')
xlim([0 100])
subplot(4,1,2)
plot(psedo_space,SFS.v,'r-')
hold on 
plot(psedo_space,playa_filt.Hz_20.v(n_start:n_end,tower_height),'k-')
plot(psedo_space,Resolved.v,'b-')
grid on
ylabel('v [ms$^{-1}$]')
xlim([0 100])
subplot(4,1,3)
plot(psedo_space,SFS.w,'r-')
hold on 
plot(psedo_space,playa_filt.Hz_20.w(n_start:n_end,tower_height),'k-')
plot(psedo_space,Resolved.w,'b-')
grid on
ylabel('w [ms$^{-1}$]')
xlim([0 100])
subplot(4,1,4)
plot(psedo_space,SFS.T(:,tower_height),'r-')
hold on 
plot(psedo_space,playa_filt.Hz_20.T(n_start:n_end,tower_height),'k-')
plot(psedo_space,Resolved.T(:,tower_height),'b-')
grid on
ylabel('T [$^{\circ}$C]')
xlim([0 100])

%% Compute u_star and tke & L 

Resolved.u_star = ((Resolved.u.*Resolved.w).^2+(Resolved.v.*Resolved.w).^2).^(1/4);
Resolved.tke = 0.5.*(Resolved.u.^2.+Resolved.v.^2+Resolved.w.^2);
Resolved.L = (Resolved.u_star.^3.*Resolved.Th)./(0.4*9.8.*Resolved.wTh);
Resolved.min_30.u_star=((Resolved.min_30.uw).^2+(Resolved.min_30.vw.^2)).^(1/4);
Resolved.min_30.tke = 0.5.*(Resolved.min_30.uu.^2.+Resolved.min_30.vv.^2+Resolved.min_30.ww.^2);
Resolved.min_30.L = (Resolved.min_30.u_star.^3.*Resolved.min_30.Th)./(0.4*9.8.*Resolved.min_30.wTh);

playa_filt.Hz_20.u_star = ((playa_filt.Hz_20.u(n_start:n_end,tower_height).*playa_filt.Hz_20.w(n_start:n_end,tower_height)).^2+...
    (playa_filt.Hz_20.v(n_start:n_end,tower_height).*playa_filt.Hz_20.w(n_start:n_end,tower_height)).^2).^(1/4);
playa_filt.Hz_20.tke = 0.5.*(playa_filt.Hz_20.u(n_start:n_end,tower_height).^2.+...
    playa_filt.Hz_20.v(n_start:n_end,tower_height).^2+playa_filt.Hz_20.w(n_start:n_end,tower_height).^2);
playa_filt.Hz_20.L = (playa_filt.Hz_20.u_star.^3.*playa_filt.Hz_20.Th(n_start:n_end,tower_height))./(0.4*9.8.*playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.Th(n_start:n_end,tower_height));
playa_filt.min_30.u_star = ((playa_filt.min_30.uw_prime.^2+playa_filt.min_30.vw_prime).^2).^(1/4);
playa_filt.min_30.tke = 0.5.*(playa_filt.min_30.uu_prime+ playa_filt.min_30.vv_prime+playa_filt.min_30.ww_prime);

%% look at TKE and u_star and L
figure()
subplot(3,1,1)
h1 = plot(psedo_space,SFS.u_star,'r-');
hold on
h2 = plot(psedo_space,Resolved.u_star,'b-');
h3 = plot(psedo_space,playa_filt.Hz_20.u_star,'k-');
legend([h1 h2 h3],'SFS','Resolved','Total')
xlim([0 20])
ylabel('$u_*$')
xlabel('x [m]')
grid on

subplot(3,1,2)
plot(psedo_space,sqrt(SFS.tke),'r-')
hold on
plot(psedo_space,sqrt(Resolved.tke),'b-')
plot(psedo_space,sqrt(playa_filt.Hz_20.tke),'k-')
xlim([0 20])
ylabel('$\sqrt{e}$')
xlabel('x [m]')
grid on

subplot(3,1,3)
plot(psedo_space,SFS.L,'r-')
hold on
plot(psedo_space,Resolved.L,'b-')
plot(psedo_space,playa_filt.Hz_20.L,'k-')
xlim([0 20])
ylabel('$L$')
xlabel('x [m]')
grid on

%% time avg variables at different times 
freq = 20; 
z = 2.02;
Pr = 0.3;
z0 = 0.11e-3; %roughness length [m] for Playa (Jensen et al)
%30 minutes
% Resolved.min_30.w = sliding_average(Resolved.w,freq,30*60);
% Resolved.min_30.v = sliding_average(Resolved.v,freq,30*60);
% Resolved.min_30.u = sliding_average(Resolved.u,freq,30*60);
% Resolved.min_30.Th = sliding_average(Resolved.Th,freq,30*60);
% Resolved.min_30.uw = sliding_average(Resolved.u.*Resolved.w,freq,30*60);
% Resolved.min_30.vw = sliding_average(Resolved.v.*Resolved.w,freq,30*60);
% Resolved.min_30.wT = sliding_average(Resolved.w.*Resolved.T(:,tower_height),freq,30*60);
% Resolved.min_30.wTh = sliding_average(Resolved.w.*Resolved.Th,freq,30*60);
% Resolved.min_30.uu = sliding_average(Resolved.u.*Resolved.u,freq,30*60);
% Resolved.min_30.vv = sliding_average(Resolved.v.*Resolved.v,freq,30*60);
% Resolved.min_30.ww = sliding_average(Resolved.w.*Resolved.w,freq,30*60);
% 
% playa_filt.min_30.uw_prime = slidingplaya_filt.Hz_20.u_prime(n_start:n_end,tower_height).*...
%     playa_filt.Hz_20.w_prime(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));
% 
% playa_filt.min_30.vw_prime = squeeze(mean(reshape(playa_filt.Hz_20.v_prime(n_start:n_end,tower_height).*...
%     playa_filt.Hz_20.v_prime(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));
% 
% playa_filt.min_30.wT = squeeze(mean(reshape(playa_filt.Hz_20.w(n_start:n_end,tower_height).*...
%     playa_filt.Hz_20.T(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));
% 
% playa_filt.min_30.wTh = squeeze(mean(reshape(playa_filt.Hz_20.w(n_start:n_end,tower_height).*...
%     playa_filt.Hz_20.Th(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));
% 
% playa_filt.min_30.uu_prime = squeeze(mean(reshape(playa_filt.Hz_20.u_prime(n_start:n_end,tower_height).*...
%     playa_filt.Hz_20.u_prime(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));
% 
% playa_filt.min_30.vv_prime = squeeze(mean(reshape(playa_filt.Hz_20.v_prime(n_start:n_end,tower_height).*...
%     playa_filt.Hz_20.v_prime(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));
% 
% playa_filt.min_30.ww_prime = squeeze(mean(reshape(playa_filt.Hz_20.w_prime(n_start:n_end,tower_height).*...
%     playa_filt.Hz_20.w_prime(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));


playa_filt.min_30.L = sliding_average(playa_filt.Hz_20.L,freq,30*60);
playa_filt.min_15.L = sliding_average(playa_filt.Hz_20.L,freq,15*60);
playa_filt.min_5.L = sliding_average(playa_filt.Hz_20.L,freq,5*60);
playa_filt.min_1.L = sliding_average(playa_filt.Hz_20.L,freq,1*60);
playa_filt.sec_30.L = sliding_average(playa_filt.Hz_20.L,freq,30);
playa_filt.sec_15.L = sliding_average(playa_filt.Hz_20.L,freq,15);
playa_filt.sec_5.L = sliding_average(playa_filt.Hz_20.L,freq,5);
playa_filt.sec_1.L = sliding_average(playa_filt.Hz_20.L,freq,1);

SFS.min_30.w = sliding_average(SFS.w,freq,30*60);
SFS.min_30.v = sliding_average(SFS.v,freq,30*60);
SFS.min_30.u = sliding_average(SFS.u,freq,30*60);
SFS.min_30.T = sliding_average(SFS.T,freq,30*60);
SFS.min_30.Th = sliding_average(SFS.Th,freq,30*60);
SFS.min_30.uw = sliding_average(SFS.u.*SFS.w,freq,30*60);
SFS.min_30.vw = sliding_average(SFS.v.*SFS.w,freq,30*60);
SFS.min_30.wT = sliding_average(SFS.w.*SFS.T(:,tower_height)',freq,30*60);
SFS.min_30.wTh = sliding_average(SFS.w.*SFS.Th,freq,30*60);
SFS.min_30.uu = sliding_average(SFS.u.*SFS.u,freq,30*60);
SFS.min_30.vv = sliding_average(SFS.v.*SFS.v,freq,30*60);
SFS.min_30.ww = sliding_average(SFS.w.*SFS.w,freq,30*60);
SFS.min_30.u_star=((SFS.min_30.uw).^2+(SFS.min_30.vw.^2)).^(1/4);
SFS.min_30.tke = 0.5.*(SFS.min_30.uu+SFS.min_30.vv+SFS.min_30.ww);
SFS.min_30.L = (SFS.min_30.u_star.^3.*SFS.min_30.Th)./(0.4*9.8.*SFS.min_30.wTh);
SFS.min_30.MOST_ra = ra_fun(z,SFS.min_30.tke,Pr,z0,'MOST',SFS.min_30.u_star,z./SFS.min_30.L); 
SFS.min_30.Shao_ra = ra_fun(z,SFS.min_30.tke,Pr,z0,'Shao',SFS.min_30.u_star,z./SFS.min_30.L);
playa_filt.min_30.T_sur = sliding_average(playa_filt.Hz_20.T_sur,freq,30*60);
playa_filt.min_30.T = sliding_average(playa_filt.Hz_20.T(n_start:n_end,tower_height),freq,30*60);
psedo_space_30min = sliding_average(psedo_space,freq,30*60);
fprintf(' ~30 min complete \n')

SFS.min_15.w = sliding_average(SFS.w,freq,60*15);
SFS.min_15.v = sliding_average(SFS.v,freq,60*15);
SFS.min_15.u = sliding_average(SFS.u,freq,60*15);
SFS.min_15.T = sliding_average(SFS.T,freq,60*15);
SFS.min_15.Th = sliding_average(SFS.Th,freq,60*15);
SFS.min_15.uw = sliding_average(SFS.u.*SFS.w,freq,60*15);
SFS.min_15.vw = sliding_average(SFS.v.*SFS.w,freq,60*15);
SFS.min_15.wT = sliding_average(SFS.w.*SFS.T(:,tower_height)',freq,60*15);
SFS.min_15.wTh = sliding_average(SFS.w.*SFS.Th,freq,60*15);
SFS.min_15.uu = sliding_average(SFS.u.*SFS.u,freq,60*15);
SFS.min_15.vv = sliding_average(SFS.v.*SFS.v,freq,60*15);
SFS.min_15.ww = sliding_average(SFS.w.*SFS.w,freq,60*15);
SFS.min_15.u_star=((SFS.min_15.uw).^2+(SFS.min_15.vw.^2)).^(1/4);
SFS.min_15.tke = 0.5.*(SFS.min_15.uu+SFS.min_15.vv+SFS.min_15.ww);
SFS.min_15.L = (SFS.min_15.u_star.^3.*SFS.min_15.Th)./(0.4*9.8.*SFS.min_15.wTh);
SFS.min_15.MOST_ra = ra_fun(z,SFS.min_15.tke,Pr,z0,'MOST',SFS.min_15.u_star,z./SFS.min_15.L); 
SFS.min_15.Shao_ra = ra_fun(z,SFS.min_15.tke,Pr,z0,'Shao',SFS.min_15.u_star,z./SFS.min_15.L);
playa_filt.min_15.T_sur = sliding_average(playa_filt.Hz_20.T_sur,freq,15*60);
playa_filt.min_15.T = sliding_average(playa_filt.Hz_20.T(n_start:n_end,tower_height),freq,15*60);
psedo_space_15min = sliding_average(psedo_space,freq,15*60);
fprintf(' ~15 min complete \n')

SFS.min_5.w = sliding_average(SFS.w,freq,60*5);
SFS.min_5.v = sliding_average(SFS.v,freq,60*5);
SFS.min_5.u = sliding_average(SFS.u,freq,60*5);
SFS.min_5.T = sliding_average(SFS.T,freq,60*5);
SFS.min_5.Th = sliding_average(SFS.Th,freq,60*5);
SFS.min_5.uw = sliding_average(SFS.u.*SFS.w,freq,60*5);
SFS.min_5.vw = sliding_average(SFS.v.*SFS.w,freq,60*5);
SFS.min_5.wT = sliding_average(SFS.w.*SFS.T(:,tower_height)',freq,60*5);
SFS.min_5.wTh = sliding_average(SFS.w.*SFS.Th,freq,60*5);
SFS.min_5.uu = sliding_average(SFS.u.*SFS.u,freq,60*5);
SFS.min_5.vv = sliding_average(SFS.v.*SFS.v,freq,60*5);
SFS.min_5.ww = sliding_average(SFS.w.*SFS.w,freq,60*5);
SFS.min_5.u_star=((SFS.min_5.uw).^2+(SFS.min_5.vw.^2)).^(1/4);
SFS.min_5.tke = 0.5.*(SFS.min_5.uu+SFS.min_5.vv+SFS.min_5.ww);
SFS.min_5.L = (SFS.min_5.u_star.^3.*SFS.min_5.Th)./(0.4*9.8.*SFS.min_5.wTh);
SFS.min_5.MOST_ra = ra_fun(z,SFS.min_5.tke,Pr,z0,'MOST',SFS.min_5.u_star,z./SFS.min_5.L); 
SFS.min_5.Shao_ra = ra_fun(z,SFS.min_5.tke,Pr,z0,'Shao',SFS.min_5.u_star,z./SFS.min_5.L);
playa_filt.min_5.T_sur = sliding_average(playa_filt.Hz_20.T_sur,freq,5*60);
playa_filt.min_5.T= sliding_average(playa_filt.Hz_20.T(n_start:n_end,tower_height),freq,5*60);
psedo_space_5min = sliding_average(psedo_space,freq,5*60);
fprintf(' ~5 min complete \n')

SFS.min_1.w = sliding_average(SFS.w,freq,60*1);
SFS.min_1.v = sliding_average(SFS.v,freq,60*1);
SFS.min_1.u = sliding_average(SFS.u,freq,60*1);
SFS.min_1.T = sliding_average(SFS.T,freq,60*1);
SFS.min_1.Th = sliding_average(SFS.Th,freq,60*1);
SFS.min_1.uw = sliding_average(SFS.u.*SFS.w,freq,60*1);
SFS.min_1.vw = sliding_average(SFS.v.*SFS.w,freq,60*1);
SFS.min_1.wT = sliding_average(SFS.w.*SFS.T(:,tower_height)',freq,60*1);
SFS.min_1.wTh = sliding_average(SFS.w.*SFS.Th,freq,60*1);
SFS.min_1.uu = sliding_average(SFS.u.*SFS.u,freq,60*1);
SFS.min_1.vv = sliding_average(SFS.v.*SFS.v,freq,60*1);
SFS.min_1.ww = sliding_average(SFS.w.*SFS.w,freq,60*1);
SFS.min_1.u_star=((SFS.min_1.uw).^2+(SFS.min_1.vw.^2)).^(1/4);
SFS.min_1.tke = 0.5.*(SFS.min_1.uu+SFS.min_1.vv+SFS.min_1.ww);
SFS.min_1.L = (SFS.min_1.u_star.^3.*SFS.min_1.Th)./(0.4*9.8.*SFS.min_1.wTh);
SFS.min_1.MOST_ra = ra_fun(z,SFS.min_1.tke,Pr,z0,'MOST',SFS.min_1.u_star,z./SFS.min_1.L); 
SFS.min_1.Shao_ra = ra_fun(z,SFS.min_1.tke,Pr,z0,'Shao',SFS.min_1.u_star,z./SFS.min_1.L);
playa_filt.min_1.T_sur = sliding_average(playa_filt.Hz_20.T_sur,freq,1*60);
playa_filt.min_1.T = sliding_average(playa_filt.Hz_20.T(n_start:n_end,tower_height),freq,1*60);
psedo_space_1min = sliding_average(psedo_space,freq,1*60);
fprintf(' ~1 min complete \n')

SFS.sec_30.w = sliding_average(SFS.w,freq,30*1);
SFS.sec_30.v = sliding_average(SFS.v,freq,30*1);
SFS.sec_30.u = sliding_average(SFS.u,freq,30*1);
SFS.sec_30.T = sliding_average(SFS.T,freq,30*1);
SFS.sec_30.Th = sliding_average(SFS.Th,freq,30*1);
SFS.sec_30.uw = sliding_average(SFS.u.*SFS.w,freq,30*1);
SFS.sec_30.vw = sliding_average(SFS.v.*SFS.w,freq,30*1);
SFS.sec_30.wT = sliding_average(SFS.w.*SFS.T(:,tower_height)',freq,30*1);
SFS.sec_30.wTh = sliding_average(SFS.w.*SFS.Th,freq,30*1);
SFS.sec_30.uu = sliding_average(SFS.u.*SFS.u,freq,30*1);
SFS.sec_30.vv = sliding_average(SFS.v.*SFS.v,freq,30*1);
SFS.sec_30.ww = sliding_average(SFS.w.*SFS.w,freq,30*1);
SFS.sec_30.u_star=((SFS.sec_30.uw).^2+(SFS.sec_30.vw.^2)).^(1/4);
SFS.sec_30.tke = 0.5.*(SFS.sec_30.uu+SFS.sec_30.vv+SFS.sec_30.ww);
SFS.sec_30.L = (SFS.sec_30.u_star.^3.*SFS.sec_30.Th)./(0.4*9.8.*SFS.sec_30.wTh);
SFS.sec_30.MOST_ra = ra_fun(z,SFS.sec_30.tke,Pr,z0,'MOST',SFS.sec_30.u_star,z./SFS.sec_30.L); 
SFS.sec_30.Shao_ra = ra_fun(z,SFS.sec_30.tke,Pr,z0,'Shao',SFS.sec_30.u_star,z./SFS.sec_30.L);
playa_filt.sec_30.T_sur = sliding_average(playa_filt.Hz_20.T_sur,freq,30*1);
playa_filt.sec_30.T = sliding_average(playa_filt.Hz_20.T(n_start:n_end,tower_height),freq,30*1);
psedo_space_30sec = sliding_average(psedo_space,freq,30*1);
fprintf(' ~30 sec complete \n')

SFS.sec_15.w = sliding_average(SFS.w,freq,15*1);
SFS.sec_15.v = sliding_average(SFS.v,freq,15*1);
SFS.sec_15.u = sliding_average(SFS.u,freq,15*1);
SFS.sec_15.T = sliding_average(SFS.T,freq,15*1);
SFS.sec_15.Th = sliding_average(SFS.Th,freq,15*1);
SFS.sec_15.uw = sliding_average(SFS.u.*SFS.w,freq,15*1);
SFS.sec_15.vw = sliding_average(SFS.v.*SFS.w,freq,15*1);
SFS.sec_15.wT = sliding_average(SFS.w.*SFS.T(:,tower_height)',freq,15*1);
SFS.sec_15.wTh = sliding_average(SFS.w.*SFS.Th,freq,15*1);
SFS.sec_15.uu = sliding_average(SFS.u.*SFS.u,freq,15*1);
SFS.sec_15.vv = sliding_average(SFS.v.*SFS.v,freq,15*1);
SFS.sec_15.ww = sliding_average(SFS.w.*SFS.w,freq,15*1);
SFS.sec_15.u_star=((SFS.sec_15.uw).^2+(SFS.sec_15.vw.^2)).^(1/4);
SFS.sec_15.tke = 0.5.*(SFS.sec_15.uu+SFS.sec_15.vv+SFS.sec_15.ww);
SFS.sec_15.L = (SFS.sec_15.u_star.^3.*SFS.sec_15.Th)./(0.4*9.8.*SFS.sec_15.wTh);
SFS.sec_15.MOST_ra = ra_fun(z,SFS.sec_15.tke,Pr,z0,'MOST',SFS.sec_15.u_star,z./SFS.sec_15.L); 
SFS.sec_15.Shao_ra = ra_fun(z,SFS.sec_15.tke,Pr,z0,'Shao',SFS.sec_15.u_star,z./SFS.sec_15.L);
playa_filt.sec_15.T_sur = sliding_average(playa_filt.Hz_20.T_sur,freq,15*1);
playa_filt.sec_15.T = sliding_average(playa_filt.Hz_20.T(n_start:n_end,tower_height),freq,15*1);
psedo_space_15sec = sliding_average(psedo_space,freq,15);
fprintf(' ~15 sec complete \n')



SFS.sec_5.w = sliding_average(SFS.w,freq,5*1);
SFS.sec_5.v = sliding_average(SFS.v,freq,5*1);
SFS.sec_5.u = sliding_average(SFS.u,freq,5*1);
SFS.sec_5.T = sliding_average(SFS.T,freq,5*1);
SFS.sec_5.Th = sliding_average(SFS.Th,freq,5*1);
SFS.sec_5.uw = sliding_average(SFS.u.*SFS.w,freq,5*1);
SFS.sec_5.vw = sliding_average(SFS.v.*SFS.w,freq,5*1);
SFS.sec_5.wT = sliding_average(SFS.w.*SFS.T(:,tower_height)',freq,5*1);
SFS.sec_5.wTh = sliding_average(SFS.w.*SFS.Th,freq,5*1);
SFS.sec_5.uu = sliding_average(SFS.u.*SFS.u,freq,5*1);
SFS.sec_5.vv = sliding_average(SFS.v.*SFS.v,freq,5*1);
SFS.sec_5.ww = sliding_average(SFS.w.*SFS.w,freq,5*1);
SFS.sec_5.u_star=((SFS.sec_5.uw).^2+(SFS.sec_5.vw.^2)).^(1/4);
SFS.sec_5.tke = 0.5.*(SFS.sec_5.uu+SFS.sec_5.vv+SFS.sec_5.ww);
SFS.sec_5.L = (SFS.sec_5.u_star.^3.*SFS.sec_5.Th)./(0.4*9.8.*SFS.sec_5.wTh);
SFS.sec_5.MOST_ra = ra_fun(z,SFS.sec_5.tke,Pr,z0,'MOST',SFS.sec_5.u_star,z./SFS.sec_5.L); 
SFS.sec_5.Shao_ra = ra_fun(z,SFS.sec_5.tke,Pr,z0,'Shao',SFS.sec_5.u_star,z./SFS.sec_5.L);
playa_filt.sec_5.T_sur = sliding_average(playa_filt.Hz_20.T_sur,freq,5*1);
playa_filt.sec_5.T = sliding_average(playa_filt.Hz_20.T(n_start:n_end,tower_height),freq,5*1);
psedo_space_5sec = sliding_average(psedo_space,freq,5);
fprintf(' ~5 sec complete \n')



SFS.sec_1.w = sliding_average(SFS.w,freq,1*1);
SFS.sec_1.v = sliding_average(SFS.v,freq,1*1);
SFS.sec_1.u = sliding_average(SFS.u,freq,1*1);
SFS.sec_1.T = sliding_average(SFS.T,freq,1*1);
SFS.sec_1.Th = sliding_average(SFS.Th,freq,1*1);
SFS.sec_1.uw = sliding_average(SFS.u.*SFS.w,freq,1*1);
SFS.sec_1.vw = sliding_average(SFS.v.*SFS.w,freq,1*1);
SFS.sec_1.wT = sliding_average(SFS.w.*SFS.T(:,tower_height)',freq,1*1);
SFS.sec_1.wTh = sliding_average(SFS.w.*SFS.Th,freq,1*1);
SFS.sec_1.uu = sliding_average(SFS.u.*SFS.u,freq,1*1);
SFS.sec_1.vv = sliding_average(SFS.v.*SFS.v,freq,1*1);
SFS.sec_1.ww = sliding_average(SFS.w.*SFS.w,freq,1*1);
SFS.sec_1.u_star=((SFS.sec_1.uw).^2+(SFS.sec_1.vw.^2)).^(1/4);
SFS.sec_1.tke = 0.5.*(SFS.sec_1.uu+SFS.sec_1.vv+SFS.sec_1.ww);
SFS.sec_1.L = (SFS.sec_1.u_star.^3.*SFS.sec_1.Th)./(0.4*9.8.*SFS.sec_1.wTh);
SFS.sec_1.MOST_ra = ra_fun(z,SFS.sec_1.tke,Pr,z0,'MOST',SFS.sec_1.u_star,z./SFS.sec_1.L); 
SFS.sec_1.Shao_ra = ra_fun(z,SFS.sec_1.tke,Pr,z0,'Shao',SFS.sec_1.u_star,z./SFS.sec_1.L);
playa_filt.sec_1.T_sur = sliding_average(playa_filt.Hz_20.T_sur,freq,1*1);
playa_filt.sec_1.T = sliding_average(playa_filt.Hz_20.T(n_start:n_end,tower_height),freq,1*1);
psedo_space_1sec = sliding_average(psedo_space,freq,1);
fprintf(' ~1 sec complete \n')

SFS.u_star=(((SFS.u.*SFS.w).^2)+((SFS.v.*SFS.w).^2)).^(1/4);
SFS.tke = 0.5.*(SFS.u.^2 + SFS.v.^2 +SFS.w.^2);
SFS.L = (SFS.u_star.^3.*SFS.Th)./(0.4*9.8.*SFS.wTh);
SFS.MOST_ra = ra_fun(z,SFS.tke,Pr,z0,'MOST',SFS.u_star,z./SFS.L'); 
SFS.Shao_ra = ra_fun(z,SFS.tke,Pr,z0,'Shao',SFS.u_star,z./SFS.L);
fprintf(' ~instantaneous complete \n')
%% Examine the different averaging times and their effect on the time series heat flux
x_start = 5000;
x_end= 5100;
figure()
subplot(2,2,1)
plot(psedo_space,SFS.Shao_ra.^-1,'-r')
hold on
plot(psedo_space,SFS.MOST_ra.^-1,'-b')
plot(psedo_space,SFS.wT'./(SFS.T(:,tower_height)-playa_filt.Hz_20.T_sur'),'-k')
grid on 
xlim([x_start x_end])
ylabel('$r_a^{-1}$') 
title('0.05 sec')

subplot(2,2,2)
plot(psedo_space_1sec,SFS.sec_1.Shao_ra.^-1,'-r')
hold on
plot(psedo_space_1sec,SFS.sec_1.MOST_ra.^-1,'-b')
plot(psedo_space_1sec,SFS.sec_1.wT./(SFS.sec_1.T-playa_filt.sec_1.T_sur),'-k')
grid on 
xlim([x_start x_end])
ylabel('$\overline{r_a}^{-1}$')
title('1 sec')

subplot(2,2,3)
plot(psedo_space_15sec,SFS.sec_15.Shao_ra.^-1,'-r')
hold on
plot(psedo_space_15sec,SFS.sec_15.MOST_ra.^-1,'-b')
plot(psedo_space_15sec,SFS.sec_15.wT./(SFS.sec_15.T-playa_filt.sec_15.T_sur),'-k')
grid on 
xlim([x_start x_end])
ylabel('$\overline{r_a}^{-1}$')
title('15 sec')

subplot(2,2,4)
plot(psedo_space_30sec,SFS.sec_30.Shao_ra.^-1,'-r')
hold on
plot(psedo_space_30sec,SFS.sec_30.MOST_ra.^-1,'-b')
plot(psedo_space_30sec,SFS.sec_30.wT./(SFS.sec_30.T-playa_filt.sec_30.T_sur),'-k')
grid on 
xlim([x_start x_end])
ylabel('$\overline{r_a}^{-1}$')
title('30 sec')





%% compute STD from "Truth"
%HERE
std_array = [.05,1,5,15,30;0,0,0,0,0;0,0,0,0,0];
%Row 2 is MOST
std_array(2,1) = sqrt((1/length(SFS.MOST_ra)).*sum((SFS.MOST_ra'.^-1 - (SFS.wT'./(playa_filt.Hz_20.T(n_start:n_end,tower_height)-playa_filt.Hz_20.T_sur'))).^2));
std_array(2,2) = sqrt((1/length(SFS.sec_1.MOST_ra)).*sum((SFS.sec_1.MOST_ra.^-1 - (SFS.sec_1.wT./(playa_filt.sec_1.T-playa_filt.sec_1.T_sur))).^2));
std_array(2,3) = sqrt((1/length(SFS.sec_5.MOST_ra)).*sum((SFS.sec_5.MOST_ra.^-1 - (SFS.sec_5.wT./(playa_filt.sec_5.T-playa_filt.sec_5.T_sur))).^2));
std_array(2,4) = sqrt((1/length(SFS.sec_15.MOST_ra)).*sum((SFS.sec_15.MOST_ra.^-1 - (SFS.sec_15.wT./(playa_filt.sec_15.T-playa_filt.sec_15.T_sur))).^2));
std_array(2,5) = sqrt((1/length(SFS.sec_30.MOST_ra)).*sum((SFS.sec_30.MOST_ra.^-1 - (SFS.sec_30.wT./(playa_filt.sec_30.T-playa_filt.sec_30.T_sur))).^2));
% std_array(2,4) = sqrt((1/length(SFS.min_1.MOST_ra)).*sum((SFS.min_1.MOST_ra.^-1 - (SFS.min_1.wT./(playa_filt.min_1.T-playa_filt.min_1.T_sur))).^2));
% std_array(2,5) = sqrt((1/length(SFS.min_5.MOST_ra)).*sum((SFS.min_5.MOST_ra.^-1 - (SFS.min_5.wT./(playa_filt.min_5.T-playa_filt.min_5.T_sur))).^2));
% std_array(2,6) = sqrt((1/length(SFS.min_15.MOST_ra)).*sum((SFS.min_15.MOST_ra.^-1 - (SFS.min_15.wT./(playa_filt.min_15.T-playa_filt.min_15.T_sur))).^2));
% std_array(2,7) = sqrt((1/length(SFS.min_30.MOST_ra)).*sum((SFS.min_30.MOST_ra.^-1 - (SFS.min_30.wT./(playa_filt.min_30.T-playa_filt.min_30.T_sur))).^2));

%Row 3 is Shao
std_array(3,1) = sqrt((1/length(SFS.Shao_ra)).*sum((SFS.Shao_ra'.^-1 - (SFS.wT'./(playa_filt.Hz_20.T(n_start:n_end,tower_height)-playa_filt.Hz_20.T_sur'))).^2));
std_array(3,2) = sqrt((1/length(SFS.sec_1.Shao_ra)).*sum((SFS.sec_1.Shao_ra.^-1 - (SFS.sec_1.wT./(playa_filt.sec_1.T-playa_filt.sec_1.T_sur))).^2));
std_array(3,3) = sqrt((1/length(SFS.sec_5.Shao_ra)).*sum((SFS.sec_5.Shao_ra.^-1 - (SFS.sec_5.wT./(playa_filt.sec_5.T-playa_filt.sec_5.T_sur))).^2));
std_array(3,4) = sqrt((1/length(SFS.sec_15.Shao_ra)).*sum((SFS.sec_15.Shao_ra.^-1 - (SFS.sec_15.wT./(playa_filt.sec_15.T-playa_filt.sec_15.T_sur))).^2));
std_array(3,5) = sqrt((1/length(SFS.sec_30.Shao_ra)).*sum((SFS.sec_30.Shao_ra.^-1 - (SFS.sec_30.wT./(playa_filt.sec_30.T-playa_filt.sec_30.T_sur))).^2));
% std_array(3,4) = sqrt((1/length(SFS.min_1.Shao_ra)).*sum((SFS.min_1.Shao_ra.^-1 - (SFS.min_1.wT./(playa_filt.min_1.T-playa_filt.min_1.T_sur))).^2));
% std_array(3,5) = sqrt((1/length(SFS.min_5.Shao_ra)).*sum((SFS.min_5.Shao_ra.^-1 - (SFS.min_5.wT./(playa_filt.min_5.T-playa_filt.min_5.T_sur))).^2));
% std_array(3,6) = sqrt((1/length(SFS.min_15.Shao_ra)).*sum((SFS.min_15.Shao_ra.^-1 - (SFS.min_15.wT./(playa_filt.min_15.T-playa_filt.min_15.T_sur))).^2));
% std_array(3,7) = sqrt((1/length(SFS.min_30.Shao_ra)).*sum((SFS.min_30.Shao_ra.^-1 - (SFS.min_30.wT./(playa_filt.min_30.T-playa_filt.min_30.T_sur))).^2));
%
figure()
h1 = loglog(std_array(1,:),std_array(2,:),'r-o');
h1.MarkerSize = 8;
hold on
h2 = loglog(std_array(1,:),std_array(3,:),'k-o');
h2.MarkerSize = 8;
grid on
legend('MOST','Shao $\textit{et al.}$')
axis tight
ylabel('$\sqrt{ \frac{1}{N} \sum(\frac{1}{r_a^{model}}- \frac{1}{r_a^{expected}})^2 }$')
xlabel('time (seconds)')

%% Compute aerodynamic resistances, (Constant surface thermal roughness and Pr number)

%Use 30 minute average (need to possibly comput on 30 min average)
zeta = z./[ones(1,(20*30*60)).*playa_filt.min_30.L(74,tower_height) ...
    ones(1,20*30*60).*playa_filt.min_30.L(75,tower_height)...
    ones(1,20*60*30).*playa_filt.min_30.L(76,tower_height) ];



SFS.ra.MOST = ra_fun(z,SFS.tke,Pr,z0,'MOST',SFS.u_star,zeta); 
SFS.ra.Shao = ra_fun(z,SFS.tke,Pr,z0,'Shao',SFS.u_star,zeta);
Resolved.ra.MOST = ra_fun(z,Resolved.tke,Pr,z0,'MOST',Resolved.u_star,zeta); 
Resolved.ra.Shao = ra_fun(z,Resolved.tke,Pr,z0,'Shao',Resolved.u_star,zeta);
playa_filt.Hz_20.ra.MOST = ra_fun(z,playa_filt.Hz_20.tke',Pr,z0,'MOST',playa_filt.Hz_20.u_star,zeta'); 
playa_filt.Hz_20.ra.Shao = ra_fun(z,playa_filt.Hz_20.tke',Pr,z0,'Shao',playa_filt.Hz_20.u_star,zeta'); 

SFS.min_30.MOST_ra = ra_fun(z,SFS.min_30.tke,Pr,z0,'MOST',SFS.min_30.u_star',(z./playa_filt.min_30.L(chunk_30min_start:chunks_30min_end-1,tower_height))); 
SFS.min_30.Shao_ra = ra_fun(z,SFS.min_30.tke,Pr,z0,'Shao',SFS.min_30.u_star,(z./playa_filt.min_30.L(chunk_30min_start:chunks_30min_end-1,tower_height))); 
Resolved.min_30.MOST_ra = ra_fun(z,Resolved.min_30.tke,Pr,z0,'MOST',Resolved.min_30.u_star',(z./playa_filt.min_30.L(chunk_30min_start:chunks_30min_end-1,tower_height))); 
Resolved.min_30.Shao_ra = ra_fun(z,Resolved.min_30.tke,Pr,z0,'Shao',Resolved.min_30.u_star,(z./playa_filt.min_30.L(chunk_30min_start:chunks_30min_end,tower_height))); 
playa_filt.min_30.ra.MOST = ra_fun(z,playa_filt.min_30.tke',Pr,z0,'MOST',playa_filt.min_30.u_star',(z./playa_filt.min_30.L(chunk_30min_start:chunks_30min_end-1,tower_height))); 
playa_filt.min_30.ra.Shao = ra_fun(z,playa_filt.min_30.tke',Pr,z0,'Shao',playa_filt.min_30.u_star,(z./playa_filt.min_30.L(chunk_30min_start:chunks_30min_end,tower_height))); 


%% Compare Ra^-1 for each method ~instantaneous
figure()
subplot(2,2,1)
h1 = plot(psedo_space,((playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end,tower_height))./...
    (playa_filt.Hz_20.T(n_start:n_end,tower_height)-playa_filt.Hz_20.T_sur')),'r-');
h1.LineWidth = 1;
hold on 
h2 = plot(psedo_space,((Resolved.w.*Resolved.T(:,tower_height)')'./(Resolved.T(:,tower_height)-playa_filt.Hz_20.T_sur')),'r-.');
h2.LineWidth = 1;
h3 = plot(psedo_space,(SFS.wT'./(SFS.T(:,tower_height)-playa_filt.Hz_20.T_sur')),'r--' );
h3.LineWidth = 1;
h4 = plot(psedo_space,playa_filt.Hz_20.ra.Shao.^(-1),'k-');
h4.LineWidth = 1;
h4.MarkerSize = 3;
h5 =plot(psedo_space,Resolved.ra.Shao.^(-1),'k-.');
h5.LineWidth = 1;
h5.MarkerSize = 3;
h6 = plot(psedo_space,SFS.ra.Shao.^(-1),'k--');
h6.LineWidth = 1;
%plot(linspace(0,100,10),linspace(0,0,10),'k-')
legend([h1 h2 h3 h4 h5 h6],'$wT/T_a-T_0$','$\widetilde{w}\widetilde{T}/\widetilde{T_a}-T_0$',...
    '$w^{\prime}T^{\prime}/T_{a}^{\prime}-T_0$','Tower','Resolved','SFS')
xlabel('x [m]')
ylabel('ra$^{-1}$ [ms$^{-1}$]')
ylim([-.6 1])
xlim([0 20])

grid on
title('Shao \it{et al.}')

subplot(2,2,2)
h1 = plot(psedo_space,((playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end,tower_height))./...
    (playa_filt.Hz_20.T(n_start:n_end,tower_height)-playa_filt.Hz_20.T_sur')),'r-');
h1.LineWidth = 1;
hold on 
h2 = plot(psedo_space,((Resolved.w.*Resolved.T(:,tower_height)')'./(Resolved.T(:,tower_height)-playa_filt.Hz_20.T_sur')),'r-.');
h2.LineWidth = 1;
h3 = plot(psedo_space,(SFS.wT'./(SFS.T(:,tower_height)-playa_filt.Hz_20.T_sur')),'r--' );
h3.LineWidth = 1;
h4 = plot(psedo_space,playa_filt.Hz_20.ra.MOST.^(-1),'k-');
h4.LineWidth = 1;
h4.MarkerSize = 3;
h5 =plot(psedo_space,Resolved.ra.MOST.^(-1),'k-.');
h5.LineWidth = 1;
h5.MarkerSize = 3;
h6 = plot(psedo_space,SFS.ra.MOST.^(-1),'k--');
h6.LineWidth = 1;
%plot(linspace(0,100,10),linspace(0,0,10),'k-')
xlabel('x [m]')
ylim([-.6 1])
xlim([0 20])
grid on
title('MOST')
% Average the ra over 30 min
mk_size = 10; 

subplot(2,2,3)
h1 = plot(psedo_space(18000:36000:end),(squeeze(mean(reshape((playa_filt.Hz_20.w(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.T(n_start:n_end,tower_height)),[36000,3]),1,'omitnan'))./...
    (playa_filt.min_30.T(chunk_30min_start:chunks_30min_end-1,5)'-playa_filt.min_30.T_sur)),'r-o');
h1.LineWidth = 1;
h1.MarkerSize = mk_size;
hold on 
h2 = plot(psedo_space(18000:36000:end),(Resolved.min_30.wT)./...
    (squeeze(mean(reshape(Resolved.T(:,tower_height),[36000,3]),1,'omitnan'))-playa_filt.min_30.T_sur),'r-.o');
h2.LineWidth = 1;
h2.MarkerSize = mk_size;
h3 = plot(psedo_space(18000:36000:end),(SFS.min_30.wT)./...
    (squeeze(mean(reshape(SFS.T(:,tower_height),[36000,3]),1,'omitnan'))-playa_filt.min_30.T_sur),'r--o');
h3.LineWidth = 1;
h3.MarkerSize = mk_size;
h4 = plot(psedo_space(18000:36000:end),playa_filt.min_30.ra.Shao.^(-1),'k-o');
h4.MarkerSize = mk_size;
h4.LineWidth = 1;
h5 =plot(psedo_space(18000:36000:end),Resolved.min_30.Shao_ra.^(-1),'k-.o');
h5.MarkerSize = mk_size;
h5.LineWidth = 1;
h6 = plot(psedo_space(18000:36000:end),SFS.min_30.Shao_ra.^(-1),'k--o');
h6.MarkerSize = mk_size;
h6.LineWidth = 1;
%plot(linspace(0,1000,10),linspace(0,0,10),'k-')
xlabel('x [m]')
ylabel('$\overline{ra}_{30min}^{-1}$ [ms$^{-1}$]')
%xlim([0 5755])
%ylim([10e-6 10^2])
axis tight
grid on
ylim([-.6 1])

subplot(2,2,4)
h1 = plot(psedo_space(18000:36000:end),(squeeze(mean(reshape((playa_filt.Hz_20.w(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.T(n_start:n_end,tower_height)),[36000,3]),1,'omitnan'))./...
    (playa_filt.min_30.T(chunk_30min_start:chunks_30min_end-1,5)'-playa_filt.min_30.T_sur)),'r-o');
h1.LineWidth = 1;
h1.MarkerSize = mk_size;
hold on 
h2 = plot(psedo_space(18000:36000:end),(Resolved.min_30.wT)./...
    (squeeze(mean(reshape(Resolved.T(:,tower_height),[36000,3]),1,'omitnan'))-playa_filt.min_30.T_sur),'r-.o');
h2.LineWidth = 1;
h2.MarkerSize = mk_size;
h3 = plot(psedo_space(18000:36000:end),(SFS.min_30.wT)./...
    (squeeze(mean(reshape(SFS.T(:,tower_height),[36000,3]),1,'omitnan'))-playa_filt.min_30.T_sur),'r--o');
h3.LineWidth = 1;
h3.MarkerSize = mk_size;
h4 = plot(psedo_space(18000:36000:end),playa_filt.min_30.ra.MOST.^(-1),'k-o');
h4.MarkerSize = mk_size;
h4.LineWidth = 1;
h5 =plot(psedo_space(18000:36000:end),Resolved.min_30.MOST_ra.^(-1),'k-.o');
h5.MarkerSize = mk_size;
h5.LineWidth = 1;
h6 = plot(psedo_space(18000:36000:end),SFS.min_30.MOST_ra.^(-1),'k--o');
h6.MarkerSize = mk_size;
h6.LineWidth = 1;
%semilogy(linspace(0,1000,10),linspace(0,0,10),'k-')
xlabel('x [m]')
%ylabel('$\overline{ra}_{30min}^{-1}$ [ms$^{-1}$]')
%xlim([0 5755])
%ylim([10e-6 10^2])
axis tight
grid on
ylim([-.6 1])
%

%% Sanity Check (wT)_sgs = (w'T')_sgs+(w)_sgs(T)_sgs

tmp = 200;
figure()
h1 = plot(psedo_space(1:tmp:end),squeeze(mean(reshape(playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end,tower_height),[tmp,180*3]),...
     1,'omitnan')),'k-');
h1.LineWidth = 3;
hold on
h3 = plot(psedo_space(1:tmp:end), squeeze(mean(reshape(SFS.wT,[tmp,180*3]),...
     1,'omitnan')),'k--');
h4 = plot(psedo_space(1:tmp:end), squeeze(mean(reshape(Resolved.w.*Resolved.T,[tmp,180*3]),...
     1,'omitnan')),'k-.');
plot(linspace(0,100,10),linspace(0,0,10),'k-')
h2=plot(psedo_space(1:tmp:end),squeeze(mean(reshape( Resolved.w.*Resolved.T+SFS.wT,[tmp,180*3]),...
     1,'omitnan')),'r--');
h2.LineWidth = 3;
legend([h1 h3 h4 h2],'$wT$','$(wT)_{sfs}$','$\widetilde{w}\widetilde{T}$','$(wT)_{sfs}+\widetilde{w}\widetilde{T}$'...
    ,'Location','north','Orientation','horizontal')
legend('boxoff')
grid on
title('2$\Delta$ Averaging')
axis([0 20 -15 15])
xlabel('x [m]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Contribution  fraction to the flux
heights = [25.5 19 10 5 2.02 0.61];
figure()
plot(percent_resolved,heights,'k-*')
hold on
plot(percent_SFS,heights,'k--o')
plot(linspace(1,1,10),linspace(0,25.5,10),'k-')
grid on
legend('$\frac{1}{N} \int_{0}^{N} \frac{\widetilde{w}\widetilde{T}}{\widetilde{wT}} dx$','$\frac{1}{N} \int_{0}^{N}  \frac{(wT)_{sfs}}{\widetilde{wT}} dx$')
axis tight
ylabel('z [m]')
yticks([0.61,2.02,5,10,19,25.5]')
%%
rho_Cp_air = 1.1614.*1.007e3; % [J kg^-1 K^-1]300k Heat Trans book
SFS.shao_flux = (SFS.T_interp(:,20)-SFS.T_interp(:,1))./SFS.ra.Shao';
SFS.MOST_flux = (SFS.T_interp(:,20)-SFS.T_interp(:,1))./SFS.ra.MOST';
Resolved.shao_flux = (Resolved.T_interp(:,20)-Resolved.T_interp(:,1))./Resolved.ra.MOST';
Resolved.MOST_flux = (Resolved.T_interp(:,20)-Resolved.T_interp(:,1))./Resolved.ra.MOST';
playa_filt.Hz_20.shao_flux = (playa_filt.Hz_20.T_interp(:,20)-playa_filt.Hz_20.T_interp(:,1))./playa_filt.Hz_20.ra.Shao';
playa_filt.Hz_20.MOST_flux = (playa_filt.Hz_20.T_interp(:,20)-playa_filt.Hz_20.T_interp(:,1))./playa_filt.Hz_20.ra.MOST;
 Sanity Check (wT)_sgs = (w'T')_sgs+(w)_sgs(T)_sgs
figure()
h1 = plot(psedo_space,rho_Cp_air.*playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end,tower_height),...
    'k-');
h1.LineWidth = 3;
hold on
h2 = plot(psedo_space,rho_Cp_air.*SFS.wT,'k--' );
h3 = plot(psedo_space, rho_Cp_air.*Resolved.w.*Resolved.T(:,tower_height)','k-.');
plot(linspace(0,100,10),linspace(0,0,10),'k-')
h4 = plot(psedo_space, rho_Cp_air.*SFS.shao_flux,'b--');
hold on
h5 = plot(psedo_space,rho_Cp_air.*SFS.MOST_flux,'r--');
h6 = plot(psedo_space,rho_Cp_air.*Resolved.shao_flux,'b-^');
h7 = plot(psedo_space,rho_Cp_air.*Resolved.MOST_flux,'r-^');
h6 = plot(psedo_space,rho_Cp_air.*playa_filt.Hz_20.shao_flux,'b-o');
h7 = plot(psedo_space,rho_Cp_air.*playa_filt.Hz_20.MOST_flux,'r-o');
h6=plot(psedo_space, rho_Cp_air.*(Resolved.w.*Resolved.T(:,tower_height)'+SFS.wT),'r--');
h6.LineWidth = 3;
legend([h1 h2 h3 h4,h5],'$\widetilde{wT}$','$(wT)_{sfs}$','$\widetilde{w}\widetilde{T}$','$H_{Shao}$','$H_{MOST}$'...
   ,'Location','north','Orientation','horizontal')
legend('boxoff')
grid on
xlim([80 100])
axis([0 20 -25 25])
xlabel('x [m]')

%% Compute Temperature gradient and Pr number 
% Get surface T from Langre Poly fit for T. 
% Potentials for heat flux are computed from resolved field!!!

dt = 1/20;
%P = zeros(501,1);
dz_2 = .2;
for i = 1:length(Resolved.T)
    %fit to poly
    Resolved.T_interp(i,:) = lagrangepoly([5,2.02,0.61],Resolved.T(i,4:6),0:.1:5);
    %Resolved.u_interp(i,:) = lagrangepoly([5,2.02,0.61],Resolved.u(i,4:6),0:.1:5);
    playa_filt.Hz_20.T_interp(i,:) = lagrangepoly([5,2.02,0.61],playa_filt.Hz_20.T(i,4:6),0:.1:5);
    SFS.T_interp(i,:) = lagrangepoly([5,2.02,0.61],SFS.T(i,4:6),0:.1:5);
    %take derivate as a central difference around 2.02 m. There for with
    %interpalation of 501 points we have a 2dz = .02
    %
    %Resolved.dT_dz(i) = (Resolved.T_interp(i,21) - Resolved.T_interp(i,19))/dz_2;
   % Resolved.du_dz(i) = (Resolved.u_interp(i,21) - Resolved.u_interp(i,19))/dz_2;
    if mod(i,10000)==0
        index = i
    end
end

Pr = (SFS.u.*SFS.w.*Resolved.dT_dz)./(SFS.w.*SFS.T(:,tower_height)'.*Resolved.du_dz);
figure()
plot(psedo_space,Resolved.dT_dz)
xlim([0 2000])
