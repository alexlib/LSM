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
emiss = .95;
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
%% Computing

%% time avg variables
Resolved.min_30.w = squeeze(mean(reshape(Resolved.w,[20*30*60,3]),1,'omitnan'));
Resolved.min_30.v = squeeze(mean(reshape(Resolved.v,[20*30*60,3]),1,'omitnan'));
Resolved.min_30.u = squeeze(mean(reshape(Resolved.u,[20*30*60,3]),1,'omitnan'));
Resolved.min_30.Th = squeeze(mean(reshape(Resolved.Th,[20*30*60,3]),1,'omitnan'));
Resolved.min_30.uw = squeeze(mean(reshape(Resolved.w.*Resolved.u(:,tower_height)',[20*30*60,3]),1,'omitnan'));
Resolved.min_30.vw = squeeze(mean(reshape(Resolved.w.*Resolved.v(:,tower_height)',[20*30*60,3]),1,'omitnan'));
Resolved.min_30.wT = squeeze(mean(reshape(Resolved.w.*Resolved.T(:,tower_height)',[20*30*60,3]),1,'omitnan'));
Resolved.min_30.wTh = squeeze(mean(reshape(Resolved.wTh',[20*30*60,3]),1,'omitnan'));
Resolved.min_30.uu = squeeze(mean(reshape(Resolved.u.*Resolved.u(:,tower_height)',[20*30*60,3]),1,'omitnan'));
Resolved.min_30.vv = squeeze(mean(reshape(Resolved.v.*Resolved.v(:,tower_height)',[20*30*60,3]),1,'omitnan'));
Resolved.min_30.ww = squeeze(mean(reshape(Resolved.w.*Resolved.w(:,tower_height)',[20*30*60,3]),1,'omitnan'));

SFS.min_30.u = squeeze(mean(reshape(SFS.u,[20*30*60,3]),1,'omitnan'));
SFS.min_30.v = squeeze(mean(reshape(SFS.v,[20*30*60,3]),1,'omitnan'));
SFS.min_30.w = squeeze(mean(reshape(SFS.w,[20*30*60,3]),1,'omitnan'));
SFS.min_30.Th = squeeze(mean(reshape(SFS.Th,[20*30*60,3]),1,'omitnan'));
SFS.min_30.uw = squeeze(mean(reshape(SFS.w.*SFS.u(:,tower_height)',[20*30*60,3]),1,'omitnan'));
SFS.min_30.vw = squeeze(mean(reshape(SFS.w.*SFS.v(:,tower_height)',[20*30*60,3]),1,'omitnan'));
SFS.min_30.wT = squeeze(mean(reshape(SFS.w.*SFS.T(:,tower_height)',[20*30*60,3]),1,'omitnan'));
SFS.min_30.wTh = squeeze(mean(reshape(SFS.wTh',[20*30*60,3]),1,'omitnan'));
SFS.min_30.uu = squeeze(mean(reshape(SFS.u.*SFS.u(:,tower_height)',[20*30*60,3]),1,'omitnan'));
SFS.min_30.vv = squeeze(mean(reshape(SFS.v.*SFS.v(:,tower_height)',[20*30*60,3]),1,'omitnan'));
SFS.min_30.ww = squeeze(mean(reshape(SFS.w.*SFS.w(:,tower_height)',[20*30*60,3]),1,'omitnan'));

playa_filt.min_30.uw_prime = squeeze(mean(reshape(playa_filt.Hz_20.u_prime(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.w_prime(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));

playa_filt.min_30.vw_prime = squeeze(mean(reshape(playa_filt.Hz_20.v_prime(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.v_prime(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));

playa_filt.min_30.wT = squeeze(mean(reshape(playa_filt.Hz_20.w(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.T(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));

playa_filt.min_30.wTh = squeeze(mean(reshape(playa_filt.Hz_20.w(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.Th(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));

playa_filt.min_30.uu_prime = squeeze(mean(reshape(playa_filt.Hz_20.u_prime(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.u_prime(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));

playa_filt.min_30.vv_prime = squeeze(mean(reshape(playa_filt.Hz_20.v_prime(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.v_prime(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));

playa_filt.min_30.ww_prime = squeeze(mean(reshape(playa_filt.Hz_20.w_prime(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.w_prime(n_start:n_end,tower_height),[20*30*60,3]),1,'omitnan'));


%% Compute u_star and tke & L 
SFS.u_star = ((SFS.u.*SFS.w).^2+(SFS.v.*SFS.w).^2).^(1/4);
SFS.tke = 0.5.*(SFS.u.^2.+SFS.v.^2+SFS.w.^2);
SFS.L = (SFS.u_star.^3.*SFS.Th)./(0.4*9.8.*SFS.wTh);
SFS.min_30.u_star=((SFS.min_30.uw).^2+(SFS.min_30.vw.^2)).^(1/4);
SFS.min_30.tke = 0.5.*(SFS.min_30.uu.^2.+SFS.min_30.vv.^2+SFS.min_30.ww.^2);
SFS.min_30.L = (SFS.min_30.u_star.^3.*SFS.min_30.Th)./(0.4*9.8.*SFS.min_30.wTh);


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

% playa_filt.Hz_20.u_star = ((playa_filt.Hz_20.u_prime(n_start:n_end,tower_height).*playa_filt.Hz_20.w_prime(n_start:n_end,tower_height)).^2+...
%     (playa_filt.Hz_20.v_prime(n_start:n_end,tower_height).*playa_filt.Hz_20.w_prime(n_start:n_end,tower_height)).^2).^(1/4);
% playa_filt.Hz_20.tke = 0.5.*(playa_filt.Hz_20.u_prime(n_start:n_end,tower_height).^2.+...
%     playa_filt.Hz_20.v_prime(n_start:n_end,tower_height).^2+playa_filt.Hz_20.w_prime(n_start:n_end,tower_height).^2);
% playa_filt.Hz_20.L = (playa_filt.Hz_20.u_star.^3.*playa_filt.Hz_20.Th(n_start:n_end,tower_height))./(0.4*9.8.*playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.Th(n_start:n_end,tower_height));
% playa_filt.min_30.u_star = ((playa_filt.min_30.uw_prime.^2+playa_filt.min_30.vw_prime).^2).^(1/4);
% playa_filt.min_30.tke = 0.5.*(playa_filt.min_30.uu_prime+ playa_filt.min_30.vv_prime+playa_filt.min_30.ww_prime);

%% look at TKE and u_star and L
figure()
subplot(3,1,1)
h1 = plot(psedo_space,SFS.u_star,'k-');
hold on
h2 = plot(psedo_space,Resolved.u_star,'r-');
h3 = plot(psedo_space,playa_filt.Hz_20.u_star,'b-');
legend([h1 h2 h3],'SFS','Resolved','Total')
xlim([0 20])
ylabel('$u_*$')
xlabel('x [m]')
grid on

subplot(3,1,2)
plot(psedo_space,sqrt(SFS.tke))
hold on
plot(psedo_space,sqrt(Resolved.tke))
plot(psedo_space,sqrt(playa_filt.Hz_20.tke))
legend('SFS','Resolved','Total')
xlim([0 20])
ylabel('$\sqrt{e}$')
xlabel('x [m]')
grid on

subplot(3,1,3)
plot(psedo_space,SFS.L)
hold on
plot(psedo_space,sqrt(Resolved.L))
plot(psedo_space,sqrt(playa_filt.Hz_20.L))
legend('SFS','Resolved','Total')
xlim([0 20])
ylabel('$L$')
xlabel('x [m]')
grid on
%% Compute Temperature gradient 
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
%%
Pr = (SFS.u.*SFS.w.*Resolved.dT_dz)./(SFS.w.*SFS.T(:,tower_height)'.*Resolved.du_dz);
figure()
plot(psedo_space,Resolved.dT_dz)
xlim([0 2000])


%%
z = 2.02;
%Use 30 minute average (need to possibly comput on 30 min average)
zeta = z./[ones(1,(20*30*60)).*playa_filt.min_30.L(74,tower_height) ...
    ones(1,20*30*60).*playa_filt.min_30.L(75,tower_height)...
    ones(1,20*60*30).*playa_filt.min_30.L(76,tower_height) ];
Pr = 0.3;
z0 = 0.11e-3; %roughness length [m] for Playa (Jensen et al)
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
%% Compare Ra for each method ~instantaneous
figure()
subplot(1,2,1)
h1 = plot(psedo_space,((playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end,tower_height))./...
    (playa_filt.Hz_20.T(n_start:n_end,tower_height)-playa_filt.Hz_20.T_sur')).^(-1),'r-');
h1.LineWidth = 1;
hold on 
h2 = plot(psedo_space,((Resolved.w.*Resolved.T(:,tower_height)')'./(Resolved.T(:,tower_height)-playa_filt.Hz_20.T_sur')).^(-1),'r-.');
h2.LineWidth = 1;
h3 = plot(psedo_space,(SFS.wT'./(SFS.T(:,tower_height)-playa_filt.Hz_20.T_sur')).^(-1),'r--' );
h3.LineWidth = 1;
h4 = plot(psedo_space,SFS.ra.Shao,'k-');
h4.LineWidth = 1;
h5 =plot(psedo_space,Resolved.ra.Shao,'k-.');
h5.LineWidth = 1;
h5.MarkerSize = 3;
h6 = plot(psedo_space,playa_filt.Hz_20.ra.Shao,'k-^');
h6.LineWidth = 1;
h6.MarkerSize = 3;
plot(linspace(0,100,10),linspace(0,0,10),'k-')
legend([h1 h2 h3 h4 h5 h6],'$T_a-T_0/\widetilde{wT}$','$\widetilde{T_a}-T_0/\widetilde{w}\widetilde{T}$',...
    '$(T_a)_{SFS}-T_0/(wT)_{SFS}$','SFS','Resolved','Tower')
xlabel('x [m]')
ylabel('ra [sm$^{-1}$]')
xlim([0 20])
ylim([-200 200])
grid on
title('Shao et al.')

subplot(1,2,2)
h1 = plot(psedo_space,((playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end,tower_height))./...
    (playa_filt.Hz_20.T(n_start:n_end,tower_height)-playa_filt.Hz_20.T_sur')).^(-1),'r-');
h1.LineWidth = 1;
hold on 
h2 = plot(psedo_space,((Resolved.w.*Resolved.T(:,tower_height)')'./(Resolved.T(:,tower_height)-playa_filt.Hz_20.T_sur')).^(-1),'r-.');
h2.LineWidth = 1;
h3 = plot(psedo_space,(SFS.wT'./(SFS.T(:,tower_height)-playa_filt.Hz_20.T_sur')).^(-1),'r--' );
h3.LineWidth = 1;
h4 = plot(psedo_space,SFS.ra.MOST,'k-');
h4.LineWidth = 1;
h5 =plot(psedo_space,Resolved.ra.MOST,'k-.');
h5.LineWidth = 1;
h5.MarkerSize = 3;
h6 = plot(psedo_space,playa_filt.Hz_20.ra.MOST,'k-^');
h6.LineWidth = 1;
h6.MarkerSize = 3;
plot(linspace(0,100,10),linspace(0,0,10),'k-')
xlabel('x [m]')
xlim([0 20])
ylim([-200 200])
grid on
title('MOST')


%% Compare Ra for each method ~instantaneous
figure()
subplot(2,2,1)
h1 = semilogy(psedo_space,((playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end,tower_height))./...
    (playa_filt.Hz_20.T(n_start:n_end,tower_height)-playa_filt.Hz_20.T_sur')),'r-');
h1.LineWidth = 1;
hold on 
h2 = semilogy(psedo_space,((Resolved.w.*Resolved.T(:,tower_height)')'./(Resolved.T(:,tower_height)-playa_filt.Hz_20.T_sur')),'r-.');
h2.LineWidth = 1;
h3 = semilogy(psedo_space,(SFS.wT'./(SFS.T(:,tower_height)-playa_filt.Hz_20.T_sur')),'r--' );
h3.LineWidth = 1;
h4 = semilogy(psedo_space,playa_filt.Hz_20.ra.Shao.^(-1),'k-');
h4.LineWidth = 1;
h4.MarkerSize = 3;
h5 =semilogy(psedo_space,Resolved.ra.Shao.^(-1),'k-.');
h5.LineWidth = 1;
h5.MarkerSize = 3;
h6 = semilogy(psedo_space,SFS.ra.Shao.^(-1),'k--');
h6.LineWidth = 1;
semilogy(linspace(0,100,10),linspace(0,0,10),'k-')
legend([h1 h2 h3 h4 h5 h6],'$\widetilde{wT}/T_a-T_0$','$\widetilde{w}\widetilde{T}/\widetilde{T_a}-T_0$',...
    '$w^{\prime}T^{\prime}/T_{a}^{\prime}-T_0$','Tower','Resolved','SFS')
xlabel('x [m]')
ylabel('ra$^{-1}$ [ms$^{-1}$]')
xlim([0 20])
ylim([10e-6 10^2])
grid on
title('Shao \it{et al.}')

subplot(2,2,2)
h1 = semilogy(psedo_space,((playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end,tower_height))./...
    (playa_filt.Hz_20.T(n_start:n_end,tower_height)-playa_filt.Hz_20.T_sur')),'r-');
h1.LineWidth = 1;
hold on 
h2 = semilogy(psedo_space,((Resolved.w.*Resolved.T(:,tower_height)')'./(Resolved.T(:,tower_height)-playa_filt.Hz_20.T_sur')),'r-.');
h2.LineWidth = 1;
h3 = semilogy(psedo_space,(SFS.wT'./(SFS.T(:,tower_height)-playa_filt.Hz_20.T_sur')),'r--' );
h3.LineWidth = 1;
h4 = semilogy(psedo_space,playa_filt.Hz_20.ra.MOST.^(-1),'k-');
h4.LineWidth = 1;
h4.MarkerSize = 3;
h5 =semilogy(psedo_space,Resolved.ra.MOST.^(-1),'k-.');
h5.LineWidth = 1;
h5.MarkerSize = 3;
h6 = semilogy(psedo_space,SFS.ra.MOST.^(-1),'k--');
h6.LineWidth = 1;
semilogy(linspace(0,100,10),linspace(0,0,10),'k-')
xlabel('x [m]')
ylim([10e-6 10^2])
xlim([0 20])
grid on
title('MOST')
% Average the ra over 30 min
mk_size = 10; 

subplot(2,2,3)
h1 = semilogy(psedo_space(1:36000:end),(squeeze(mean(reshape((playa_filt.Hz_20.w(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.T(n_start:n_end,tower_height)),[36000,3]),1,'omitnan'))./...
    (playa_filt.min_30.T(chunk_30min_start:chunks_30min_end-1,5)'-playa_filt.min_30.T_sur)),'r-');
h1.LineWidth = 1;
h1.MarkerSize = mk_size;
hold on 
h2 = semilogy(psedo_space(1:36000:end),(squeeze(mean(reshape((Resolved.w(:,tower_height).*...
    Resolved.T(:,tower_height)),[36000,3]),1,'omitnan'))./...
    (squeeze(mean(reshape(Resolved.T_interp(:,20),[36000,3]),1,'omitnan'))-playa_filt.min_30.T_sur)),'r-.');
h2.LineWidth = 1;
h2.MarkerSize = mk_size;
h3 = semilogy(psedo_space(1:36000:end),(squeeze(mean(reshape((SFS.w(:,tower_height).*...
    SFS.T(:,tower_height)),[36000,3]),1,'omitnan'))./...
    (squeeze(mean(reshape(SFS.T_interp(:,20),[36000,3]),1,'omitnan'))-playa_filt.min_30.T_sur)),'r--');
h3.LineWidth = 1;
h3.MarkerSize = mk_size;
h4 = semilogy(psedo_space(1:36000:end),playa_filt.min_30.ra.Shao.^(-1),'k-');
h4.MarkerSize = mk_size;
h4.LineWidth = 1;
h5 =semilogy(psedo_space(1:36000:end),Resolved.min_30.Shao_ra.^(-1),'k-.');
h5.MarkerSize = mk_size;
h5.LineWidth = 1;
h6 = semilogy(psedo_space(1:36000:end),SFS.min_30.Shao_ra.^(-1),'k--');
h6.MarkerSize = mk_size;
h6.LineWidth = 1;
semilogy(linspace(0,1000,10),linspace(0,0,10),'k-')
xlabel('x [m]')
ylabel('$\overline{ra}_{30min}^{-1}$ [ms$^{-1}$]')
xlim([0 5755])
ylim([10e-6 10^2])
grid on


subplot(2,2,4)
h1 = semilogy(psedo_space(1:36000:end),(squeeze(mean(reshape((playa_filt.Hz_20.w(n_start:n_end,tower_height).*...
    playa_filt.Hz_20.T(n_start:n_end,tower_height)),[36000,3]),1,'omitnan'))./...
    (playa_filt.min_30.T(chunk_30min_start:chunks_30min_end-1,5)'-playa_filt.min_30.T_sur)),'r-');
h1.LineWidth = 1;
h1.MarkerSize = mk_size;
hold on 
h2 = semilogy(psedo_space(1:36000:end),(squeeze(mean(reshape((Resolved.w(:,tower_height).*...
    Resolved.T(:,tower_height)),[36000,3]),1,'omitnan'))./...
    (squeeze(mean(reshape(Resolved.T_interp(:,20),[36000,3]),1,'omitnan'))-playa_filt.min_30.T_sur)),'r-.');
h2.LineWidth = 1;
h2.MarkerSize = mk_size;
h3 = semilogy(psedo_space(1:36000:end),(squeeze(mean(reshape((SFS.w(:,tower_height).*...
    SFS.T(:,tower_height)),[36000,3]),1,'omitnan'))./...
    (squeeze(mean(reshape(SFS.T_interp(:,20),[36000,3]),1,'omitnan'))-playa_filt.min_30.T_sur)),'r--');
h3.LineWidth = 1;
h3.MarkerSize = mk_size;
h4 = semilogy(psedo_space(1:36000:end),playa_filt.min_30.ra.MOST.^(-1),'k-');
h4.MarkerSize = mk_size;
h4.LineWidth = 1;
h5 =semilogy(psedo_space(1:36000:end),Resolved.min_30.MOST_ra.^(-1),'k-.');
h5.MarkerSize = mk_size;
h5.LineWidth = 1;
h6 = semilogy(psedo_space(1:36000:end),SFS.min_30.MOST_ra.^(-1),'k--');
h6.MarkerSize = mk_size;
h6.LineWidth = 1;
semilogy(linspace(0,1000,10),linspace(0,0,10),'k-')
% legend([h1 h2 h3 h4 h5 h6],'$T_a-T_0/\widetilde{wT}$','$\widetilde{T_a}-T_0/\widetilde{w}\widetilde{T}$',...
%     '$(T_a)_{SFS}-T_0/(wT)_{SFS}$','Tower','Resolved','SFS')
xlabel('x [m]')
xlim([0 5755])
ylim([10e-6 10^2])
grid on
%


%%
rho_Cp_air = 1.1614.*1.007e3; % [J kg^-1 K^-1]300k Heat Trans book
SFS.shao_flux = (SFS.T_interp(:,20)-SFS.T_interp(:,1))./SFS.ra.Shao';
SFS.MOST_flux = (SFS.T_interp(:,20)-SFS.T_interp(:,1))./SFS.ra.MOST';
Resolved.shao_flux = (Resolved.T_interp(:,20)-Resolved.T_interp(:,1))./Resolved.ra.MOST';
Resolved.MOST_flux = (Resolved.T_interp(:,20)-Resolved.T_interp(:,1))./Resolved.ra.MOST';
playa_filt.Hz_20.shao_flux = (playa_filt.Hz_20.T_interp(:,20)-playa_filt.Hz_20.T_interp(:,1))./playa_filt.Hz_20.ra.Shao';
playa_filt.Hz_20.MOST_flux = (playa_filt.Hz_20.T_interp(:,20)-playa_filt.Hz_20.T_interp(:,1))./playa_filt.Hz_20.ra.MOST;
%% Sanity Check (wT)_sgs = (w'T')_sgs+(w)_sgs(T)_sgs
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
figure()
subplot(4,1,1)
plot(psedo_space,playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end, tower_height),'k')
ylabel('$wT$')
grid on
axis tight

subplot(4,1,2)
plot(psedo_space,Resolved.w.*Resolved.w,'k')
ylabel('$\tilde{w} \tilde{T}$')
grid on
axis tight

subplot(4,1,3)
plot(psedo_space,SFS.w.*SFS.w,'k')
ylabel('$w^\prime T^\prime$')
grid on 
axis tight

subplot(4,1,4)
plot(psedo_space,(SFS.w.*Resolved.T)+(Resolved.w.*SFS.T),'k')
ylabel('$w^\prime \tilde{T} + \tilde{w} T^\prime$')
xlabel('x [m]')
grid on
axis tight
%% avg 
tmp = 800;
figure()
subplot(4,1,1)
plot(psedo_space(1:tmp:end),squeeze(mean(reshape(playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end, tower_height),[tmp,45*3]),...
    1,'omitnan')),'k')
ylabel('$wT$')
grid on
axis tight

subplot(4,1,2)
plot(psedo_space(1:tmp:end),squeeze(mean(reshape(Resolved.w.*Resolved.w,[tmp,45*3]),...
    1,'omitnan')),'k')
ylabel('$\tilde{w} \tilde{T}$')
grid on
axis tight

subplot(4,1,3)
plot(psedo_space(1:tmp:end),squeeze(mean(reshape(SFS.w.*SFS.w,[tmp,45*3]),...
    1,'omitnan')),'k')
ylabel('$w^\prime T^\prime$')
grid on 
axis tight

subplot(4,1,4)
plot(psedo_space(1:tmp:end),squeeze(mean(reshape((SFS.w.*Resolved.T)+(Resolved.w.*SFS.T),[tmp,45*3]),...
    1,'omitnan')),'k')
ylabel('$w^\prime \tilde{T} + \tilde{w} T^\prime$')
xlabel('x [m]')
grid on
axis tight
%%
figure()
 plot(psedo_space(1:tmp:end),squeeze(mean(reshape(playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end, tower_height),[tmp,45*3]),...
     1,'omitnan')),'k')
 hold on
plot(psedo_space(1:tmp:end),squeeze(mean(reshape(Resolved.wT,[tmp,45*3]),...
     1,'omitnan')),'b')
plot(psedo_space(1:tmp:end),squeeze(mean(reshape(Resolved.wT-Resolved.w.*Resolved.T,[tmp,45*3]),...
     1,'omitnan')),'g')
plot(psedo_space(1:tmp:end),squeeze(mean(reshape(SFS.wT,[tmp,45*3]),...
     1,'omitnan')),'r')
%%
tmp = 800;
figure()
% plot(psedo_space(1:tmp:end),squeeze(mean(reshape(playa_filt.Hz_20.w(n_start:n_end,tower_height).*playa_filt.Hz_20.T(n_start:n_end, tower_height),[tmp,45*3]),...
%     1,'omitnan')),'k')
plot(psedo_space(1:tmp:end),squeeze(mean(reshape(Resolved.wT,[tmp,45*3]),...
     1,'omitnan')),'k')
hold on
plot(psedo_space(1:tmp:end),squeeze(mean(reshape(Resolved.w.*Resolved.w,[tmp,45*3]),...
    1,'omitnan')),'b')
plot(psedo_space(1:tmp:end),squeeze(mean(reshape(SFS.w.*SFS.w,[tmp,45*3]),...
    1,'omitnan')),'g')
plot(psedo_space(1:tmp:end),squeeze(mean(reshape((SFS.w.*Resolved.T)+(Resolved.w.*SFS.T),[tmp,45*3]),...
    1,'omitnan')),'r')
legend('$wT$','$\tilde{w} \tilde{T}$','$w^\prime T^\prime$','$w^\prime \tilde{T} + \tilde{w} T^\prime$')
grid on
axis tight
xlabel('x [m]')
