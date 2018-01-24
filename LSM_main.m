% LSM model
% GOAL: drive a coupled LSM with experimental data from MATERHORN and
% Cabow and compare the results with the Shao formualtion for the SHF to
% the MOST formualtion against experimental data

clear; clc;% close all;

%SIMULATION OPTIONS

%pick the data set to drive the simulation: ie 'MATERHORN' or 'Cabow'.
driving_data_set= 'MATERHORN';
%SHF formualtion: ie 'MOST' or 'Shao'
HF_option = 'MOST';
%plots on?: opt, 'on' or 'off'
plots = 'on';
%atmosphere respond?
atmrespondTF = 'on';
%boundary layer respondonse
ABLTF ='on';
%Days to simulate (days)
t_day = 3;
%tower height index 1 for 25m...6 for 0.5 m
tower_height = 2;


ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);

%% Constants
%Variables and constants
Cp=1005.7; Cv=719; %heat capacities @ constant pressure & volume [J/kg/K] (Appendix 2 of Emanuel [1994])
g=9.80665; %standard surface gravity [m/s2]
Rd=287.04; %Ideal Gas Constant of DRY air [J/kg/K] (Appendix 2 of Emanuel [1994])
Rv=461.40; %Ideal Gas Constant of water vapor [J/kg/K] (Appendix A.1.4 of Jacobson [1999])
sigma = 5.67e-8; %[W/m^2K^4] stephan boltzmann constant
rho_air = 1; %[kg/m^3]

%time
dt= 60;           %model timestep [s]
tmax= t_day*24*3600;  %maximum time [s]
t_hr = 0:24;

% Load data
switch(driving_data_set)
    case 'MATERHORN'
        %tower heights
       % z = [0.61 2.02 5 10.4 19 25.5]; %m
        z = [25.5 19 10.4 5 2.02 0.61];
        
        %load radiation data: 5 min data
        data_start = 6696; %25 May 0600 UTC
        data_end = data_start + (12*24); %26 May 0800 UTC
        load('./Materhorn_data/MATERHORN_Rad_data.mat');
        SWdn = interp1(linspace(0,24,289),rad_data(data_start:data_end,6),linspace(0,24,1440),'spline');
        SWup = interp1(linspace(0,24,289),rad_data(data_start:data_end,7),linspace(0,24,1440),'spline');
        LWdn = interp1(linspace(0,24,289),rad_data(data_start:data_end,8),linspace(0,24,1440),'spline');
        LWup = interp1(linspace(0,24,289),rad_data(data_start:data_end,9),linspace(0,24,1440),'spline');
        albedo = mean(SWup./SWdn,'omitnan');
        %Plot data to sneak a peak
        if isequal(plots,'on')
            figure()
            plot(linspace(0,24,1440),SWdn)
            hold on
            plot(linspace(0,24,1440),LWdn)
            plot(linspace(0,24,1440),SWup)
            legend('SWin','LWin','SWup')
            title('Driving radation data')
            ylabel('E Wm$^{-2}$')
            axis tight
        end
        % clear rad_data;
        
        %Soil property inputs
        emiss = 0.93; %Jin & Liang [2006]
        z0 = 0.11e-3; %roughness length [m] for Playa (Jensen et al)
        T0 = 290.15;   %Deep soil temperature [K] for playa (Morrison et al. 2017)
        D = 0.06; %Thermally active layer from exp data from playa (Morrison et al)
        Cp_soil = 1542.8;  %specific heat of soil organic material [J/kg/K]
        rho_soil = 1425; %density of soil organic material [kg/m3]
        alpha_soil = 0.4E-6; %Thermal diffusivity for playa [m2/s]
        k_soil = 0.9;             %themal conductivity for playa [W/mK]
        Cs = Cp_soil*rho_soil*D; %heat capacity of organic soil [J/K/m2]
        
        
        %load tower data: 30 min data
        load('/Users/travismorrison/Documents/Local_Data/MATERHORN/data/tower_data/Playa_tower_raw/PlayaSpring_raw_2013_05_24.mat')
        load('./Materhorn_data/playaSpring30minLinDetUTESpac3.mat')
        start_index = 432000; %May 24 0600 UTC or 0000MST
        end_index = start_index +(20*60*60*24); %May 25 0800 UTC or 0200MST
        u_prime = rearrangeHeights(rawFlux.uPrime(start_index:end_index-1,:));
        v_prime = rearrangeHeights(rawFlux.vPrime(start_index:end_index-1,:));
        w_prime = rearrangeHeights(rawFlux.wPrime(start_index:end_index-1,:));
        chunk = 20*60*10;
        
        %compute terms on 30 min chunk
        u_prime_sq = squeeze(mean(reshape((u_prime(:,tower_height).^2)',[chunk,144]),1,'omitnan'));
        v_prime_sq = squeeze(mean(reshape((v_prime(:,tower_height).^2)',[chunk,144]),1,'omitnan'));
        w_prime_sq = squeeze(mean(reshape((w_prime(:,tower_height).^2)',[chunk,144]),1,'omitnan'));
        
        uw_prime = squeeze(mean(reshape((u_prime(:,tower_height).*w_prime(:,tower_height))',[chunk,144]),1,'omitnan'));
        vw_prime = squeeze(mean(reshape((v_prime(:,tower_height).*w_prime(:,tower_height))',[chunk,144]),1,'omitnan'));

        tke = 0.5.*(u_prime_sq+v_prime_sq+...
            w_prime_sq);
        
        
        u_star = (uw_prime.^2+vw_prime.^2).^(1/4);
        %30 min data indexes
        data_start = 1072;%24 May 0800 UTC
        data_end = data_start + (2*24);%25 May 0800 UTC
        U = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.spdAndDir(data_start:data_end,3:3:18)),linspace(0,24,1440),'spline');
        u = detrend(interp1(linspace(0,24,49),rearrangeHeights(playaSpring.rotatedSonic(data_start:data_end,1:3:18)),linspace(0,24,1440),'spline'));
        v = detrend(interp1(linspace(0,24,49),rearrangeHeights(playaSpring.rotatedSonic(data_start:data_end,2:3:18)),linspace(0,24,1440),'spline'));
        w = detrend(interp1(linspace(0,24,49),rearrangeHeights(playaSpring.rotatedSonic(data_start:data_end,3:3:18)),linspace(0,24,1440),'spline'));
        
        load('./Materhorn_data/u_star_sgs.mat')
        load('./Materhorn_data/tke_sgs.mat')
        u_star = interp1(linspace(0,24,144),u_star_sgs,linspace(0,24,1440));
        tke = interp1(linspace(0,24,144),tke_sgs,linspace(0,24,1440));
        L = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.L(data_start:data_end,2:end)),linspace(0,24,1440));
        T_air_tower = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.derivedT(data_start:data_end,2:4:22)),linspace(0,24,1440),'spline');
        %Computing drag coefficient for non-neutral conditions (per hour)
        %zeta = [0.405967591987653,	1.79075784723176,	0.0209220759296613,	1.70652034119992,	3.54795806836081,	15.9133087096185,	-3.02207506547757,	-0.341759085608449,	-1.20368229625590,	-1.28284369011724,	-1.11396080428715,	-1.53205435649209,	-1.03876262278506,	-0.905703575982586,	-0.954347088271912,	-0.661575432802616,	-0.709072759058575,	-0.502739096089809,	-0.710748095123254,	-0.368831557068724,	-0.281621327394183,	-0.116221268903314,	-0.0533745368836851,	-0.0621922059401636,	0.609550392995834];
        %zeta = interp1(linspace(0,24,25),zeta,linspace(0,24,1440),'spline');
        zeta = z./L;
        %zeta above calculated based on L values from the experimental data (less manual calculation steps in matlab)
        %Plot data to sneak a peak
        if isequal(plots,'on')
            figure()
            subplot(5,1,1)
            plot(linspace(0,24,1440),U)
            axis tight
            ylabel('U [ms$^{-1}$]')
            title('Driving tower data')
            axis tight
            subplot(5,1,2)
            plot(linspace(0,24,1440),tke)
            ylabel('e [m$^2$s$^{-2}$]')
            axis tight
            subplot(5,1,3)
            plot(linspace(0,24,1440),L)
            ylabel('L [m]')
            ylim([-10 10])
            subplot(5,1,4)
            plot(linspace(0,24,1440),u_star)
            axis tight
            ylabel('$u^*$ ms$^{-1}$')
            subplot(5,1,5)
            plot(linspace(0,24,1440),T_air_tower)
            ylabel('T [C]')
            xlabel('time [hrs]')
            axis tight
        end
        %initial vars
        H = playaSpring.H(data_start,2)*playaSpring.H(data_start,3)*playaSpring.H(data_start,15:12:75);
        LH= playaSpring.LHflux(data_start,6);
        BR = H./LH;
    case 'Cabow'
        %spell cabow right and add data set!
end
%%
%atmospheric conditions
RH=0.44; %mean 25m RH need to check time period
e=RH*satvap(squeeze(mean(T_air_tower(1,tower_height))),mean(H(:,tower_height)),BR(1))/100; %vapor pressure [hPa]
Psurf=1000;     %surface pressure [hPa]
qair=(Rd/Rv)*e/Psurf;  %specific humidity [g/g]
hmin=100;       %minimum height of atmospheric boundary layer [m]
gvmax=1/8000; %max vegetation conductance [m/s]; reciprocal of vegetation resistance
thetavM0=(T_air_tower(1,tower_height)+273.15)*(1+0.61*qair); %initial virtual potential temperature [K]; Eq. 1.5.1b of Stull [1988]
Beta=.2;       %closure hypothesis:  fraction of surface virtual potential temperature flux that determines entrainment heat flux
gamma=5/1000;   %5/1000   #slope of thetav above growing ABL [K/m]
qabove=qair/5; %specific humidity of air above ABL [g/g] changed from 5 to 1

%Initialize surface T based on upwave Longwave rad
T_init = (LWup(1)/(sigma*emiss))^(1/4);
T = T_init;
Ta = T_air_tower(1,tower_height);
z = z(tower_height);
W = 0; %subsidence rate
thetaM = Ta+273.15;
qM=qair;   %initialize with prescribed specific humidity [g/g]
thetavM=thetaM*(1+0.61*qM);   %virtual potential temperature [K];  Eq. 1.5.1b of Stull [1988]
%%
t_curr = 0;
h = hmin; %inital height of bl [m]
cnt = 0;
input_cnt = 1;
while (t_curr<tmax)
    %inciment counter
    cnt = cnt+1;
    
    %Compute Radiation budget
    LWup_curr = emiss*sigma*T^4; %[W/m^2]
    LWdn_curr = LWdn(input_cnt);
    SWdn_curr = SWdn(input_cnt);
    SWup_curr = SWup(input_cnt);
    Rnet = SWdn_curr - SWup_curr + LWdn_curr - LWup_curr;
    
    if isequal(atmrespondTF,'on')
        Ta = thetaM;   %air temperature [K] is the one from zero-order jump model
        qa = qM;
        %else
        %Ta = approx(x=as.numeric(names(Ta.c))*3600,y=Ta.c,xout=tcurr%%(24*3600))$y+273.15  #use prescribed value
        %qa = qair
    end
    
    %Compute SHF
    tke_curr = tke(input_cnt);
    ra = ra_fun(z,tke_curr,z0,HF_option,u_star(input_cnt),zeta(input_cnt,tower_height));
    H = ((Cp*rho_air)/ra)*(T-Ta);
    
    %Compute LHF
    beta=1;   %vegetation control (Jarvis type model)
    lambda=1000*latentheat(T-273.15); %latent heat of vaporization [J/kg]
    esat=satvap(T-273.15)/100;  %saturation specific humidity [hPa]
    e=qa*Psurf/(Rd/Rv); %vapor pressure [hPa]
    VPD=100*(esat-e);            %vapor pressure deficit [Pa]
    qstar=(Rd/Rv)*esat/Psurf;   %saturation specific humidity [g/g]
    rv=1/gvmax;
    LH = (lambda*rho_air/(ra+rv))*(qstar-qa);%latentheat(H,BR(1));
    
    %Compute GHF
    G = Rnet - LH - H;
    
    %update surface temperature
    dT=(G/Cs)*dt;
    T=T+dT;
    
    %if want atmosphere to respond
    %Based on "zero-order jump" or "slab" model of convective boundary layer, described in Pg. 151~155 of Garratt [1992]
    if isequal(atmrespondTF,'on')
        %update ABL-averaged q
        %calculate surface virtual heat flux
        lambda=latentheat(T);  %latent heat of vaporization [J/g]
        E=LH/lambda;   %surface moisture flux [g/m^2/s]
        F0theta=H/Cp;  %potential heat flux [K-kg/m^2/s]
        F0thetav=F0theta+0.073*lambda*E/Cp; %virtual heat flux [K-kg/m^2/s]
        if isequal(ABLTF,'on')
            Fhthetav=-1*Beta*F0thetav; %closure hypothesis (Eq. 6.15 of Garratt [1992])
            %calculate ABL growth rate
            dh.dt=(1+2*Beta)*F0thetav/(gamma*h);
            if F0thetav<=0.00
                dh.dt=0;   %ABL collapses
            end
        else
            dh.dt=0;
            Fhthetav=0;
        end
        %calculate entrainment flux of humidity
        deltaq=qabove - qM;
        Fhq=-1*rho_air*deltaq*(dh.dt-W)*1000; %entrainment flux of humidity [g/m2/s] NOTE:  assume CONSTANT air density!
        dq.dt=(E-Fhq)/(rho_air*1000*h); %change of humidity in ABL [1/s]; NOTE:  assume CONSTANT air density!
        qM=qM+dq.dt*dt;                %updated qM [g/g]
        %update ABL-averaged thetav
        dthetavM.dt=(F0thetav-Fhthetav)/h;  %change of thetav in ABL [K-kg/m^3/s]
        dthetavM.dt=dthetavM.dt/rho_air;        %[K-kg/m^3/s]=>[K/s]
        thetavM=thetavM+dthetavM.dt*dt;
        thetaM=thetavM/(1+0.61*qM);   %potential temperature, from virtual pot temp [K];  Eq. 1.5.1b of Stull [1988]
        %update ABL height
        h=h+dh.dt*dt;
        if F0thetav<=0.00 && isequal(ABLTF,'on')
            h=hmin; %override value:  ABL collapses
        end
    end
    
    input_cnt = input_cnt+1;
    if input_cnt>1440
        input_cnt =1;
    end
    t_curr = t_curr+dt;
    %save vars!
    result.time(cnt) = t_curr;
    result.SWdn(cnt) = SWdn_curr;
    result.SWup(cnt) = SWup_curr;
    result.LWup(cnt) = LWup_curr;
    result.LWdn(cnt) = LWdn_curr;
    result.Rnet(cnt) = Rnet;
    result.G(cnt) = G;
    result.H(cnt) = H;
    result.LH(cnt) = LH;
    result.ABL_h(cnt) = h;
    result.T(cnt) = T;
    result.Ta(cnt) = Ta;
end
%%
if isequal(plots,'on')
    figure()
    subplot(2,2,1)
    title('Simulated Temperature')
    plot(linspace(0,t_day,cnt),result.T,'-k')
    hold on
    plot(linspace(0,t_day,cnt),result.Ta,'--k')
    legend('$T_s$','$T_a$')
    xlabel('time [days]')
    ylabel('T [K]')
    
    subplot(2,2,2)
    title('Simulated SHF and LHF')
    plot(linspace(0,t_day,cnt),result.H,'-r')
    hold on
    plot(linspace(0,t_day,cnt),result.LH,'-b')
    legend('H','LH')
    xlabel('time [days]')
    ylabel('E [Wm$^{-2}$]')
    
    subplot(2,2,3)
    title('Simulated ABL')
    plot(linspace(0,t_day,cnt),result.ABL_h,'--k')
    hold on
    plot(linspace(0,t_day,cnt),linspace(2500,2500,cnt),'-k')
    legend('Simulated','Max Obs')
    xlabel('time [days]')
    ylabel('h [m]')
    axis([0 t_day 0 3000])
    
    subplot(2,2,4)
    title('Simulated CV budget (+ for in)')
    plot(linspace(0,t_day,cnt),result.Rnet,'-k')
    hold on
    plot(linspace(0,t_day,cnt),-result.G,'-g')
    plot(linspace(0,t_day,cnt),-result.H,'-r')
    plot(linspace(0,t_day,cnt),-result.LH,'-b')
    legend('R$_{net}$','G','H','$H_L$')
    xlabel('time [days]')
    ylabel('E [Wm$^{-2}$]')
    axis tight
    
    H_index = [15 27 75 39 51 63];
    figure()
    subplot(1,2,1)
    plot(linspace(0,24,cnt/t_day), result.H((cnt-cnt/t_day+1):cnt))
    hold on
    plot(linspace(0,24.5,49),playaSpring.H(data_start:data_end,2).*...
        playaSpring.H(data_start:data_end,3).*playaSpring.H(data_start:data_end,H_index(tower_height)))
    %15 for 25.5 m, 27 for 19.4 m, 63 for 0.61m, 75 for 10.4 m,39 for 5 m,
    %51 for 2 m
    legend('modeled','obs (25.5 m)')
    axis tight
    ylabel('H [Wm$^{-2}$]')
    xlabel('time [hrs]')
    
    subplot(1,2,2)
    plot(linspace(0,24,cnt/t_day), result.LH((cnt-cnt/t_day+1):cnt))
    hold on
    plot(linspace(0,24.5,49),playaSpring.LHflux(data_start:data_end,6))
    legend('modeled','obs (10.4 m)')
    axis tight
    ylabel('$H_L$ [Wm$^{-2}$]')
    xlabel('time [hrs]')
end




