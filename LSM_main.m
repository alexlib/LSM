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
        plot(SWin)
        %clear rad_data
        
        %Soil property inputs
        
        %load tower data: 30 min data
        data_start = 1116;%25 May 0000 UTC
        data_end = data_start + (2*24);%26 May 0000 UTC
        load('.\Materhorn_data\playaSpring30minLinDetUTESpac3.mat')
        %U = 
        %tke = 
        %L =  
    case 'Cabow'
        %spell cabow right and add data set!
end