% LSM model
% GOAL: drive a coupled LSM with experimental data from MATERHORN and
% Cabow and compare the results with the Shao formualtion for the SHF to
% the MOST formualtion against experimental data

clear; clc; close all;

%SIMULATION OPTIONS

%pick the data set to drive the simulation: ie 'MATERHORN' or 'Cabow'.
driving_data_set= 'MATERHORN';
%SHF formualtion: ie 'MOST' or 'Shao'
SHF_method = 'MOST';
%plots on?: opt, 'on' or 'off'
plots = 'on';
%Days to simulate (days)
t_day = 20; 

%% Load data
switch(driving_data_set)
    case 'MATERHORN'
        %tower heights
        z = [0.61 2.02 5 10.4 19 25.5]; %m
        
        %load radiation data: 5 min data
        data_start = 6695+24; %24 May 0800 UTC
        data_end = data_start + (12*24); %25 May 0800 UTC
        load('.\Materhorn_data\MATERHORN_Rad_data.mat');
        SWin = rad_data(data_start:data_end,6);
        SWout = rad_data(data_start:data_end,7);
        LWin = rad_data(data_start:data_end,8);
        LWout = rad_data(data_start:data_end,9);
        albedo = SWout./SWin;
        %Plot data to sneak a peak
        if isequal(plots,'on')
            figure()
            plot(linspace(0,24,289),SWin)
            hold on
            plot(linspace(0,24,289),LWin)
            legend('SWin','LWin')
            title('Driving radation data')
            ylabel('E Wm^{-2}')
            axis tight
        end
        clear rad_data
        
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
        U = rearrangeHeights(playaSpring.spdAndDir(data_start:data_end,3:3:18));
        tke = rearrangeHeights(playaSpring.tke(data_start:data_end,2:end));
        L = rearrangeHeights(playaSpring.L(data_start:data_end,2:end));
        T_air = rearrangeHeights(playaSpring.derivedT(data_start:data_end,2:4:22));
        %Plot data to sneak a peak
        if isequal(plots,'on')
            figure()
            subplot(4,1,1)
            plot(linspace(0,24,49),U)
            axis tight
            ylabel('U [ms^{-1}]')
            title('Driving tower data')
            axis tight
            subplot(4,1,2)
            plot(linspace(0,24,49),tke)
            ylabel('TKE [m^2s^{-2}]')
            axis tight
            subplot(4,1,3)
            plot(linspace(0,24,49),L)
            ylabel('L [m]')
            axis tight
            subplot(4,1,4)
            plot(linspace(0,24,49),T_air)
            ylabel('T [C]')
            xlabel('time [hrs]')
            axis tight
        end
    case 'Cabow'
        %spell cabow right and add data set!
end
%% Model parameters

%time
dt= 60;           %model timestep [s]
tmax= t_day*24*3600;  %maximum time [s]
t_hr = 0:24; 
T_cur = 0;

%Results dumper
results.T =0;
%% time loop
cnt = 0; 
while (tcurr<tmax)
   cnt = cnt+1; 
    
    
    
end
