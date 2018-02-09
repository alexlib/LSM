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

%Function loads data from the tower and filters for Northerly winds and
%Taylors hypothesis (u'/U),0.3
% [U_north_filt,u_filt,v_filt,w_filt,u_prime_filt,v_prime_filt,...
%     w_prime_filt,tke_filt,L_filt,T_prime_filt, Temp_filt, t_filt] = Preparedata();



%% Use taylors frozen theory hypothesis


freq = 20;
npoints = size(playa_filt.min_30.U,1);
chunk = freq*60*30; %in frames
time = chunk/freq; % [sec]
tower_height = 5;
les_grid_size = 10; %[m]

t = (playa_filt.min_30.U(:,tower_height)./les_grid_size).^(-1);

%
% figure()
% plot(datetime(datevec(playa_filt.min_30.t)),t','k*')
% hold on
% ylabel('$\frac{dx_{LES}}{U} = t_{filt} [seconds] $')
% xlabel('time [hrs]')
% grid on 

filter_size = t.*freq; %[npoints for flux data]
taylors_window = chunk; %[npoints = to length of application of taylors hypothesis]

% Filter data into Resolved and SGS
turbo_var_length = 36000*size(playa_filt.Hz_20.u,2);
playa_filt.Hz_20.u_prime = reshape(playa_filt.Hz_20.u_prime, [turbo_var_length,6]);
playa_filt.Hz_20.v_prime = reshape(playa_filt.Hz_20.v_prime, [turbo_var_length,6]);
playa_filt.Hz_20.w_prime = reshape(playa_filt.Hz_20.w_prime, [turbo_var_length,6]);
playa_filt.Hz_20.T_prime = reshape(playa_filt.Hz_20.T_prime, [turbo_var_length,6]);
playa_filt.Hz_20.T = reshape(playa_filt.Hz_20.T, [turbo_var_length,6]);
playa_filt.Hz_20.u = reshape(playa_filt.Hz_20.u, [turbo_var_length,6]);
playa_filt.Hz_20.v = reshape(playa_filt.Hz_20.v, [turbo_var_length,6]);
playa_filt.Hz_20.w = reshape(playa_filt.Hz_20.w, [turbo_var_length,6]);
%compute U on 20 Hz
for k = 1:6
    playa_filt.Hz_20.U = sqrt(playa_filt.Hz_20.u(:,k).^2+playa_filt.Hz_20.v(:,k).^2+playa_filt.Hz_20.w(:,k).^2);
end
%%
tmp =  turbo_var_length-sum(floor(filter_size));

LES.resolved.u_prime = zeros(1,tmp);
LES.resolved.v_prime = zeros(1,tmp);
LES.resolved.w_prime = zeros(1,tmp);
LES.resolved.T_prime= zeros(1,tmp);
LES.resolved.u = zeros(tmp,6);
LES.resolved.v = zeros(tmp,6);
LES.resolved.w = zeros(tmp,6);
LES.resolved.T = zeros(tmp,6);
LES.resolved.U = zeros(tmp,6);
%%

start_index = 1;
nchunks = turbo_var_length./36000;
for i = 1:nchunks
    end_index = floor(start_index + filter_size(i));
    for j = 1:taylors_window-filter_size(i)
        
        %since we need theta for the finite difference in height
        
        %primes
        LES.resolved.u_prime(start_index) = squeeze(mean(playa_filt.Hz_20.u_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        LES.resolved.v_prime(start_index) = squeeze(mean(playa_filt.Hz_20.v_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        LES.resolved.w_prime(start_index) = squeeze(mean(playa_filt.Hz_20.w_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
        LES.resolved.T_prime(start_index) = squeeze(mean(playa_filt.Hz_20.T_prime(start_index:end_index,tower_height),'omitnan')); %this solves for the grid resolved
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

for k = 1:6
    start_index = 1;
    for i = 1:nchunks
        end_index = floor(start_index + filter_size(i));
        for j = 1:taylors_window-filter_size(i)
            LES.resolved.T(start_index,k) = squeeze(mean(playa_filt.Hz_20.T(start_index:end_index,k),1,'omitnan')); %this solves for the grid resolved
            LES.resolved.u(start_index,k) = squeeze(mean(playa_filt.Hz_20.u(start_index:end_index,k),1,'omitnan')); %this solves for the grid resolved
            LES.resolved.v(start_index,k) = squeeze(mean(playa_filt.Hz_20.v(start_index:end_index,k),1,'omitnan')); %this solves for the grid resolved
            LES.resolved.w(start_index,k) = squeeze(mean(playa_filt.Hz_20.w(start_index:end_index,k),1,'omitnan')); %this solves for the grid resolved
%             if k == 1
%                 LES.resolved.U(start_index,k) = squeeze(mean(playa_filt.Hz_20.U(start_index:end_index,k),1,'omitnan')); %this solves for the grid resolved
%             end
            start_index = start_index+1;
            end_index = end_index+1;
        end
         if mod(i,100)==0
          chunk_num = i
        end
    end
    height = k
end
%% Redistribute points over which 
LES.resolved.u_prime = interp1(linspace(0,24,tmp),LES.resolved.u_prime,linspace(0,24,turbo_var_length));
LES.resolved.v_prime = interp1(linspace(0,24,tmp),LES.resolved.v_prime,linspace(0,24,turbo_var_length));
LES.resolved.w_prime = interp1(linspace(0,24,tmp),LES.resolved.w_prime,linspace(0,24,turbo_var_length));
LES.resolved.T_prime = interp1(linspace(0,24,tmp),LES.resolved.T_prime,linspace(0,24,turbo_var_length));
%%
for k = 1:6
     test1(:,k) = interp1(linspace(0,24,tmp),LES.resolved.T(:,k),linspace(0,24,turbo_var_length));
   test2(:,k) = interp1(linspace(0,24,tmp),LES.resolved.u(:,k),linspace(0,24,turbo_var_length));
   test3(:,k) = interp1(linspace(0,24,tmp),LES.resolved.v(:,k),linspace(0,24,turbo_var_length));
   test4(:,k) = interp1(linspace(0,24,tmp),LES.resolved.w(:,k),linspace(0,24,turbo_var_length));
end
%%
LES.resolved.T = test1;
LES.resolved.u = test2;
LES.resolved.v = test3;
LES.resolved.w = test4;
%%
for k = 1:6
    LES.resolved.U(:,k) = sqrt(LES.resolved.u(:,k).^2+LES.resolved.v(:,k).^2+LES.resolved.w(:,k).^2);
     LES.sgs.U(:,k) = sqrt(LES.sgs.u(:,k).^2+LES.sgs.v(:,k).^2+LES.sgs.w(:,k).^2);
end


%  LES.resolved.T(:,k) = interp1(linspace(0,24,tmp),LES.resolved.T(:,k),linspace(0,24,turbo_var_length));
%     LES.resolved.u(:,k) = interp1(linspace(0,24,tmp),LES.resolved.u(:,k),linspace(0,24,turbo_var_length));
%     LES.resolved.v(:,k) = interp1(linspace(0,24,tmp),LES.resolved.v(:,k),linspace(0,24,turbo_var_length));
%     LES.resolved.w(:,k) = interp1(linspace(0,24,tmp),LES.resolved.w(:,k),linspace(0,24,turbo_var_length));
%%
LES.sgs.u_prime = playa_filt.Hz_20.u_prime(:,tower_height)' - LES.resolved.u_prime;
LES.sgs.v_prime = playa_filt.Hz_20.v_prime(:,tower_height)' - LES.resolved.v_prime;
LES.sgs.w_prime = playa_filt.Hz_20.w_prime(:,tower_height)' - LES.resolved.w_prime;
LES.sgs.T_prime = playa_filt.Hz_20.T_prime(:,tower_height)' - LES.resolved.T_prime;
LES.sgs.u = playa_filt.Hz_20.u - LES.resolved.u;
LES.sgs.v = playa_filt.Hz_20.v - LES.resolved.v;
LES.sgs.w = playa_filt.Hz_20.w - LES.resolved.w;
LES.sgs.T = playa_filt.Hz_20.T - LES.resolved.T;
LES.sgs.U = playa_filt.Hz_20.U - LES.resolved.U;


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
