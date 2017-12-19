% LSM model
% GOAL: drive a coupled LSM with experimental data from MATERHORN and
% Cabow and compare the results with the Shao formualtion for the SHF to
% the MOST formualtion against experimental data

clear; clc; close all;

%SIMULATION OPTIONS

%pick the data set to drive the simulation: ie 'MATERHORN' or 'Cabow'.
driving_data_set= 'MATERHORN';
%SHF formualtion: ie 'MOST' or 'Shao'
HF_option = 'Shao';
%plots on?: opt, 'on' or 'off'
plots = 'on';
%Days to simulate (days)
t_day = 20;
%% Constants
%Variables and constants
Cp=1005.7; Cv=719; %heat capacities @ constant pressure & volume [J/kg/K] (Appendix 2 of Emanuel [1994])
g=9.80665; %standard surface gravity [m/s2]
Rd=287.04; %Ideal Gas Constant of DRY air [J/kg/K] (Appendix 2 of Emanuel [1994])
Rv=461.40; %Ideal Gas Constant of water vapor [J/kg/K] (Appendix A.1.4 of Jacobson [1999])
sigma = 5.67e-8; %[W/m^2K^4] stephan boltzmann constant
rho_air = 1; %[kg/m^3]
results.T =0; %Results dumper

%time
dt= 60;           %model timestep [s]
tmax= t_day*24*3600;  %maximum time [s]
t_hr = 0:24;
t_curr = 0;

% Load data

switch(driving_data_set)
    case 'MATERHORN'
        %tower heights
        z = [0.61 2.02 5 10.4 19 25.5]; %m
        
        %load radiation data: 5 min data
        data_start = 6695+24; %24 May 0800 UTC
        data_end = data_start + (12*24); %25 May 0800 UTC
        load('.\Materhorn_data\MATERHORN_Rad_data.mat');
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
            ylabel('E Wm^{-2}')
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
        data_start = 1072;%25 May 0000 UTC
        data_end = data_start + (2*24);%26 May 0000 UTC
        load('.\Materhorn_data\playaSpring30minLinDetUTESpac3.mat')
        U = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.spdAndDir(data_start:data_end,3:3:18)),linspace(0,24,1440),'spline');
        tke = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.tke(data_start:data_end,2:end)),linspace(0,24,1440),'spline');
        L = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.L(data_start:data_end,2:end)),linspace(0,24,1440),'spline');
        T_air_tower = interp1(linspace(0,24,49),rearrangeHeights(playaSpring.derivedT(data_start:data_end,2:4:22)),linspace(0,24,1440),'spline');
        %Plot data to sneak a peak
        if isequal(plots,'on')
            figure()
            subplot(4,1,1)
            plot(linspace(0,24,1440),U)
            axis tight
            ylabel('U [ms^{-1}]')
            title('Driving tower data')
            axis tight
            subplot(4,1,2)
            plot(linspace(0,24,1440),tke)
            ylabel('TKE [m^2s^{-2}]')
            axis tight
            subplot(4,1,3)
            plot(linspace(0,24,1440),L)
            ylabel('L [m]')
            axis tight
            subplot(4,1,4)
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
%  time loop

RH=0.44; %mean 25m RH need to check time period
e=RH*satvap(squeeze(mean(T_air_tower(1,1))),mean(H(:,1)),BR(1))/100; %vapor pressure [hPa]
Psurf=1000;     %surface pressure [hPa]
qair=(Rd/Rv)*e/Psurf;  %specific humidity [g/g]

%atmospheric conditions
hmin=100;       %minimum height of atmospheric boundary layer [m]
hmax=20000;      %max height of ABL [m]
thetavM0=(T_air_tower(1,1)+273.15)*(1+0.61*qair); %initial virtual potential temperature [K]; Eq. 1.5.1b of Stull [1988]
Beta=.2;       %closure hypothesis:  fraction of surface virtual potential temperature flux that determines entrainment heat flux
gamma=5/1000;   %5/1000   #slope of thetav above growing ABL [K/m]
qabove=qair/5; %specific humidity of air above ABL [g/g] changed from 5 to 1

T_init = (LWup(1)/(sigma*emiss))^(1/4);
T = T_init;
Ta = T_air_tower(1,1);
z = z(end);

thetaM = Ta+273.15;
qM=qair;   %initialize with prescribed specific humidity [g/g]
thetavM=thetaM*(1+0.61*qM);   %virtual potential temperature [K];  Eq. 1.5.1b of Stull [1988]
%%
cnt = 0;
input_cnt = 1;
while (tcurr<tmax)
    %inciment counter
    cnt = cnt+1;
    
    %Compute Radiation budget
    LWup_curr = emiss*sigma*T^4; %[W/m^2]
    LWdn_curr = LWdn(input_cnt);
    SWdn_curr = SWdw(input_cnt);
    SWup_curr = albdeo*SWdw(input_cnt);
    Rnet = SWdw_curr - SWup_curr + LWdw_curr - LWup_curr;
    
    if isequal(atmrespondTF,'on')
        Ta = thetaM;   %air temperature [K] is the one from zero-order jump model
        qa = qM;
        %else
        %Ta = approx(x=as.numeric(names(Ta.c))*3600,y=Ta.c,xout=tcurr%%(24*3600))$y+273.15  #use prescribed value
        %qa = qair
    end
    
    %Compute SHF
    tke_curr = tke(input_cnt);
    ra = ra_fun(z,tke,z0,HF_option);
    H = ((Cp*rho_air)/ra)*(T-T_air_tower);
    
    %Compute LHF
    LH = latentheat(H,BR(1));
    
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
        lambda=latentheat(H,BR(1));  %latent heat of vaporization [J/g]
        E=LE/lambda;   %surface moisture flux [g/m^2/s]
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
        Fhq=-1*rho*deltaq*(dh.dt-W)*1000; %entrainment flux of humidity [g/m2/s] NOTE:  assume CONSTANT air density!
        dq.dt=(E-Fhq)/(rho*1000*h); %change of humidity in ABL [1/s]; NOTE:  assume CONSTANT air density!
        qM=qM+dq.dt*dt;                %updated qM [g/g]
        %update ABL-averaged thetav
        dthetavM.dt=(F0thetav-Fhthetav)/h;  %change of thetav in ABL [K-kg/m^3/s]
        dthetavM.dt=dthetavM.dt/rho;        %[K-kg/m^3/s]=>[K/s]
        thetavM=thetavM+dthetavM.dt*dt;
        thetaM=thetavM/(1+0.61*qM);   %potential temperature, from virtual pot temp [K];  Eq. 1.5.1b of Stull [1988]
        %update ABL height
        h=h+dh.dt*dt;
        if F0thetav<=0.00 && isequal(ABLTF,'on')
            h=hmin; %override value:  ABL collapses
        end
    end
    
    input_cnt = intputcnt+1;
    if input_cnt>1440
        input_cnt =1;
    end
    t_curr = t_curr+dt;
end
