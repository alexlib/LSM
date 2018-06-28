%Loading data meeting MOST requirments
clear; close all;
% settings header
addpath('/Users/travismorrison/Documents/Code/functions_library')
ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);
clc;
%May 9th missing all heights for fw except 25m
%May 10th missing all heights for fw except 25m
%May 17th missing all heights for fw 
%may 14th and 15th work well
%day_array = ["May9","May10","May14","May15","May17","May20","May31"];
day_array = ["May9"];
avg_chunk = 20*60*30;
%
for day = 1:length(day_array)
    clearvars -except day day_array Shao_all Most_all wT_all avg_chunk;
    day_of_interest = day_array(day);
    fprintf(day_of_interest,': \n');
    %load data
    %Import data which has been properly orientated by UTESPAC (LPF on 30 minute chunks and GPF), filtered for
    %Northerly flow and meets taylors hypothesis requirment
    load('/Users/travismorrison/Documents/Local_Data/MATERHORN/data/tower_data/Tower_filtered/May_2013_Playa_Filtered_N_Taylors.mat');
    fprintf('Filtered Data Loaded: Filtered for Northerly Flow & u/U<0.3 \n')
    %
    % reshape data into time series from 30min chunks
    z = [25.5,19.4,10,5,2.02,0.61];
    
    
    playa_filt.Hz_20.w = reshape(playa_filt.Hz_20.w, [36000*580,6]);
    playa_filt.Hz_20.u_prime = reshape(playa_filt.Hz_20.u_prime, [36000*580,6]);
    playa_filt.Hz_20.v_prime = reshape(playa_filt.Hz_20.v_prime, [36000*580,6]);
    playa_filt.min_30.fwT = squeeze(mean(playa_filt.Hz_20.fwT,1));
    playa_filt.Hz_20.fwT = reshape(playa_filt.Hz_20.fwT, [36000*580,6]);
    playa_filt.Hz_20.fwTh_prime = reshape(playa_filt.Hz_20.fwTh_prime,[36000*580,6]);
    playa_filt.Hz_20.u = reshape(playa_filt.Hz_20.u, [36000*580,6]);
    playa_filt.Hz_20.v = reshape(playa_filt.Hz_20.v, [36000*580,6]);
    playa_filt.Hz_20.w_prime = reshape(playa_filt.Hz_20.w_prime, [36000*580,6]);
    
    time_vec = datetime(datevec(playa_filt.min_30.t));
    
    %
    % load radiation data: 5 min data for estimate of surface temperature
    load('./Materhorn_data/MATERHORN_Rad_data.mat');
    % From looking at data set, the entire day of May 10th was a solid day...
    % Lets check the Obukhov Length for validity of MOST
    %enter a day which has a continous chunk of afternoon data
    
    tower_height = 5; %2-m height index
    switch day_of_interest
        case 'May9' %May 9th 1600-1830 UTC
            start_i = 10;
            end_i = 15;
            rad_index =[2495:2525];
        case 'May10' %May 10th 1600-2400 UTC
            start_i = 43;
            end_i = 58;
            rad_index =[2495:2525];
        case 'May14'%May 14th 1600-2030UTC
            start_i = 93;
            end_i = 102;
            rad_index = [3935:3989];
        case 'May15'%May 15th 1600-1800 UTC
            start_i = 128;
            end_i = 132;
            rad_index = [4223:4247];
        case 'May17'%May 17th 1600-2330 UTC
            start_i = 180;
            end_i = 209;
            rad_index =[4799:4889];
        case 'May20'%Mary 20th 1600-1900 UTC
            start_i = 260;
            end_i = 266;
            rad_index =[5663:5699];
        case 'May31'%Mary 31th 1600-1830 UTC
            start_i = 488;
            end_i = 493;
            rad_index =[8831:8861];
    end
     %index for avg data
    start_iavg = start_i*avg_chunk;
    end_iavg = end_i*avg_chunk;
    
    %index for raw data
    start_i20hz = start_i*20*30*60;
    end_i20hz = end_i*20*30*60-1;
    
    %number of avg chunks
    n_chunks = floor((end_i20hz-start_i20hz+1)./avg_chunk);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load radition data
    LWup = interp1(linspace(0,24,size(rad_index,2)),rad_data(rad_index,9),linspace(0,24,(end_i20hz-start_i20hz+1)),'spline');
    SWdn = interp1(linspace(0,24,size(rad_index,2)),rad_data(rad_index,6),linspace(0,24,(end_i20hz-start_i20hz+1)),'spline');
    emiss = .97;
    sb= 5.67E-8;
    playa_filt.Hz_20.T_sur = (LWup./(emiss*sb)).^(1/4)-273.15;
    playa_filt.avg.T_sur = squeeze(mean(reshape(playa_filt.Hz_20.T_sur,[avg_chunk,n_chunks]),1,'omitnan'));
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINE AVG CHUNK VARS
    %Compute fluxes and u_star
       %constant flux layer?
    playa_filt.avg.wTh_prime = squeeze(mean(reshape(playa_filt.Hz_20.w_prime(start_i20hz:end_i20hz,:).*...
        playa_filt.Hz_20.fwTh_prime(start_i20hz:end_i20hz,:),[avg_chunk,n_chunks,6]),1,'omitnan'));
    playa_filt.avg.uw_prime = squeeze(mean(reshape(playa_filt.Hz_20.u_prime(start_i20hz:end_i20hz,:).*...
        playa_filt.Hz_20.w_prime(start_i20hz:end_i20hz,:),[avg_chunk,n_chunks,6]),1,'omitnan'));
    
    
    %Ustar and sqrt(tke)?
    for h = 1:length(z)
        playa_filt.avg.u_star(:,h) = ((squeeze(mean(reshape(playa_filt.Hz_20.u_prime(start_i20hz:end_i20hz,h).*...
            playa_filt.Hz_20.w_prime(start_i20hz:end_i20hz,h), [avg_chunk,(n_chunks)]),1,'omitnan'))).^2).^(1/4);
        playa_filt.avg.U(:,h) = squeeze(mean(reshape(playa_filt.Hz_20.u(start_i20hz:end_i20hz,h), [avg_chunk,(n_chunks)]),1,'omitnan'));
        playa_filt.avg.u_prime2(:,h) = squeeze(mean(reshape(playa_filt.Hz_20.u_prime(start_i20hz:end_i20hz,h).^2, [avg_chunk,(n_chunks)]),1,'omitnan'));
        playa_filt.avg.v_prime2(:,h) = squeeze(mean(reshape(playa_filt.Hz_20.v_prime(start_i20hz:end_i20hz,h).^2, [avg_chunk,(n_chunks)]),1,'omitnan'));
        playa_filt.avg.w_prime2(:,h) = squeeze(mean(reshape(playa_filt.Hz_20.w_prime(start_i20hz:end_i20hz,h).^2, [avg_chunk,(n_chunks)]),1,'omitnan'));
        playa_filt.avg.tke(:,h) = 0.5.*(playa_filt.avg.u_prime2(:,h)+playa_filt.avg.v_prime2(:,h)+playa_filt.avg.w_prime2(:,h));
        playa_filt.avg.fwTh(:,h)= squeeze(mean(reshape(playa_filt.Hz_20.fwT(start_i20hz:end_i20hz,h), [avg_chunk,(n_chunks)]),1,'omitnan'));
        playa_filt.avg.L(:,h) = -(playa_filt.avg.u_star(:,h).^3.*playa_filt.avg.fwTh(:,h))./(9.8.*0.4*playa_filt.avg.wTh_prime(:,h));
    end
    
    playa_filt.avg.t = squeeze(mean(reshape(playa_filt.Hz_20.t(start_i20hz:end_i20hz), [avg_chunk,(n_chunks)]),1,'omitnan'));
   
    
  
   %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure (1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %look at L
    figure('rend','painters','pos',[10 10 1200 1200])
    subplot(3,3,1)
    plot(datetime(datevec(playa_filt.avg.t)),z(tower_height)./...
        playa_filt.avg.L(:,tower_height),'k-o')
    grid on
    ylabel('$\zeta$')
    title(day_of_interest)
    
    %Subsdiance?
    subplot(3,3,2)
    plot(linspace(0,(end_i-start_i)*30,(end_i20hz-start_i20hz+1)),playa_filt.Hz_20.w(start_i20hz:end_i20hz,tower_height)','k-')
    grid on
    ylabel('w')
    xlabel('time (minutes)')
    axis tight
    
 
    
    subplot(3,3,3)
    plot(playa_filt.avg.wTh_prime,z,'k-o')
    xlabel('$\overline{w^\prime \theta^\prime}$')
    ylabel('z (m)')
    grid on
    
    subplot(3,3,4)
    plot(playa_filt.avg.uw_prime,z,'k-o')
    xlabel('$\overline{u^\prime w^\prime}$')
    ylabel('z (m)')
    grid on
    
    subplot(3,3,5)
    plot(datetime(datevec(playa_filt.avg.t)),playa_filt.avg.T_sur,'k')
    grid on
    axis tight
    ylabel('$T_s$')
   
    subplot(3,3,6)
    plot(datetime(datevec(playa_filt.Hz_20.t(start_i20hz:end_i20hz))),SWdn,'k')
    grid on
    axis tight
    ylabel('$S_{\downarrow}$')
    
    subplot(3,3,7)
    plot(datetime(datevec(playa_filt.Hz_20.t(start_i20hz:end_i20hz))),playa_filt.Hz_20.fwT(start_i20hz:end_i20hz,tower_height),'k')
    hold on
    plot(datetime(datevec(playa_filt.min_30.t(start_i:end_i))),playa_filt.min_30.T(start_i:end_i,tower_height),'r')
    grid on
    axis tight
    ylabel('$\theta_a$')
    
    subplot(3,3,8)
    plot(datetime(datevec(playa_filt.Hz_20.t(start_i20hz:end_i20hz))),playa_filt.Hz_20.fwTh_prime(start_i20hz:end_i20hz,tower_height),'k')
    hold on
    plot(datetime(datevec(playa_filt.Hz_20.t(start_i20hz:end_i20hz))),playa_filt.Hz_20.w_prime(start_i20hz:end_i20hz,tower_height),'r')
    grid on
    axis tight
    legend('$\theta^\prime$','$w^\prime$')
    

    subplot(3,3,9)
    plot(datetime(datevec(playa_filt.avg.t)),(.4.*playa_filt.avg.u_star(:,tower_height))','k-o')
    hold on
    plot(datetime(datevec(playa_filt.avg.t)),(.5.*sqrt(playa_filt.avg.tke(:,tower_height))),'r-o')
    grid on
    legend('$\kappa u_*$','$C_k \sqrt{e}/Pr$')
    axis tight
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Compute MRD on 1Hr chunks for wT & uw
    M = floor(log2(20*60*60));
    n_chunks_mrd = floor((end_i-start_i)./2);
    mrd_start_i = start_i20hz;
    mrd_end_i = start_i20hz+(20*60*60);
    
    for  h = 1:length(z)
        for n = 1:n_chunks_mrd
            uw_MRD(:,h,n) = compute_MRD(playa_filt.Hz_20.u_prime(mrd_start_i:mrd_end_i,h),...
                playa_filt.Hz_20.w_prime(mrd_start_i:mrd_end_i,h),M,0);
            wT_MRD(:,h,n) = compute_MRD(playa_filt.Hz_20.w_prime(mrd_start_i:mrd_end_i,h),...
                playa_filt.Hz_20.fwTh_prime(mrd_start_i:mrd_end_i,h),M,0);
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
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute gradients
    % %look at the T profile first over the period examined
    
    for i = 1:n_chunks
        [ fitzT, fitT(:,i), dTdz(:,i)] = nat_log_x_fit(z,playa_filt.avg.fwTh(i,:));
        [ fitzU, fitU(:,i), dUdz(:,i)] = nat_log_x_fit(z,playa_filt.avg.U(i,:));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure (3) Ln fit for T and
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure()
    subplot(1,2,1)
    plot(playa_filt.avg.fwTh,z,'k*');
    hold on
    plot(fitT,fitzT,'b-')
    axis tight
    %legend(h1, h2, 'data','ln fit')
    xlabel('T ($^{o}$C)')
    ylabel('z (m)')
    grid on
    
    subplot(1,2,2)
    plot(playa_filt.avg.U,z,'k*');
    hold on
    %plot(playa_filt.min_30.T_sur,0,'k*')
    plot(fitU,fitzU,'b-')
    axis tight
    xlabel('U (m/s)')
    grid on
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the gradient Richardson Number
    g = 9.8; %Gravitational Constant
    Rg = ((g./playa_filt.avg.fwTh(:,tower_height)).*dTdz(tower_height,:)')./dUdz(tower_height,:)'.^2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Examine various ways to look at Prt computed at 2-m
    
    Prt_case = 'direct_compute';
    switch Prt_case
        case 'constant'
            Prt = 0.3; %Turbulent Pr number from Shao et al
        case 'direct_compute'
            Prt = (playa_filt.avg.uw_prime(:,tower_height).*dTdz(tower_height,:)')./(playa_filt.avg.wTh_prime(:,tower_height).*dUdz(tower_height,:)');
        case 'Spectral'
            [ Prt_neutral, Prt ] = Compute_Prt_Spectral_approach( Rg, z(tower_height)./playa_filt.avg.L(:,tower_height) , z(tower_height));
    end
    
 %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute MOST and Shao flux on 30 minute averages
    z0 = 0.11e-3; %roughness length [m] for Playa (Jensen et al)
    
    %
    shift = 20*30*60;
    for i = 1:n_chunks
        start = start_i20hz+((i-1)*shift);
        finish = start +shift;
        [ tke_sgs,u_star_sgs, tke_resolved, u_star_resolved, q_sgs, q_resolved ] = compute_sgs_tke( playa_filt.avg.U(i,5), playa_filt.Hz_20.u(start:finish),...
            playa_filt.Hz_20.v(start:finish) ,playa_filt.Hz_20.w(start:finish),playa_filt.Hz_20.fwT(start:finish));
        
        tke_mean_sgs(i) = mean(tke_sgs,'omitnan');
        u_star_mean_sgs(i) = mean(u_star_sgs,'omitnan');
        tke_mean_resolved(i) = mean(tke_resolved,'omitnan');
        u_star_mean_resolved(i) = mean(u_star_resolved,'omitnan');
        
        q_mean_sgs(i) = mean(q_sgs,'omitnan');
        q_mean_resolved(i) = mean(q_resolved,'omitnan');
    end
    
    
    %%
  
    flux_start_i = start_i20hz;
    flux_end_i = start_i20hz+(20*60*30);
    rho_Cp_air = 1.1614.*1.007e3;
    % for  h = 1:length(z)
    %     MOST_ra(:,h) = ra_fun(z(h),playa_filt.min_30.tke(start_i:end_i-1,h),Prt,z0,'MOST',...
    %         playa_filt.min_30.u_star(:,h),z(h)./playa_filt.min_30.L(start_i:end_i-1,h));
    %     Shao_ra(:,h) = ra_fun(z(h),playa_filt.min_30.tke(start_i:end_i-1,h),Prt,z0,'Shao',...
    %         playa_filt.min_30.u_star(:,h),z(h)./playa_filt.min_30.L(start_i:end_i,h));
    %     MOST_H(:,h) = -rho_Cp_air.*((playa_filt.min_30.T(start_i:end_i-1,h) - playa_filt.min_30.T_sur')./MOST_ra(:,h));
    %     Shao_H(:,h) = -rho_Cp_air.*((playa_filt.min_30.T(start_i:end_i-1,h) - playa_filt.min_30.T_sur')./Shao_ra(:,h));
    % end
    %
    
    for  h = 1:length(z)
        %SGS
        MOST_ra(:,h) = ra_fun(z(h),tke_mean_sgs,Prt,z0,'MOST',...
            u_star_mean_sgs',z(h)./playa_filt.avg.L(:,h));
        Shao_ra(:,h) = ra_fun(z(h),tke_mean_sgs',Prt,z0,'Shao',...
            u_star_mean_sgs',z(h)./playa_filt.avg.L(:,h));
        MOST_H(:,h) = -rho_Cp_air.*((playa_filt.avg.fwTh(:,h) - playa_filt.avg.T_sur')./MOST_ra(:,h));
        Shao_H(:,h) = -rho_Cp_air.*((playa_filt.avg.fwTh(:,h) - playa_filt.avg.T_sur')./Shao_ra(:,h));
        
        %Resolved
        %      MOST_ra_resolved(:,h) = ra_fun(z(h),tke_mean_resolved,Prt,z0,'MOST',...
        %         u_star_mean_resolved',z(h)./playa_filt.min_30.L(start_i:end_i-1,h));
        %     Shao_ra_resolved(:,h) = ra_fun(z(h),tke_mean_resolved',Prt,z0,'Shao',...
        %         u_star_mean_resolved',z(h)./playa_filt.min_30.L(start_i:end_i,h));
        %     MOST_H_resolved(:,h) = -rho_Cp_air.*((playa_filt.min_30.T(start_i:end_i-1,h) - playa_filt.min_30.T_sur')./MOST_ra_resolved(:,h));
        %     Shao_H_resolved(:,h) = -rho_Cp_air.*((playa_filt.min_30.T(start_i:end_i-1,h) - playa_filt.min_30.T_sur')./Shao_ra_resolved(:,h));
    end
    
    
    
    playa_filt.avg.H = rho_Cp_air.*playa_filt.avg.wTh_prime;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure (4) Fluxes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure()
    errorbar(datenum(playa_filt.avg.t),mean(MOST_H,2,'omitnan'),...
        std(MOST_H,0,2,'omitnan'),'r-o')
    hold on
    errorbar(datenum(playa_filt.avg.t),mean(Shao_H,2,'omitnan'),...
        std(Shao_H,0,2,'omitnan'),'b-o')
    errorbar(datenum(playa_filt.avg.t),...
        mean(playa_filt.avg.H,2,'omitnan'),std(playa_filt.avg.H,0,2,'omitnan'),'k-o')
    
    errorbar(datenum(playa_filt.avg.t),...
        mean(rho_Cp_air.*q_mean_resolved',2,'omitnan'),std(rho_Cp_air.*q_mean_resolved',0,2,'omitnan'),'g-o')
    
    errorbar(datenum(playa_filt.avg.t),...
        mean(rho_Cp_air.*q_mean_sgs',2,'omitnan'),std(rho_Cp_air.*q_mean_sgs',0,2,'omitnan'),'m-o')
    
    
    legend('MOST','Shao','Tower')
    axis tight
    datetick('x','mmmdd HH:MM','keeplimits')
    
    ylabel('H (Wm$^{-2}$)')
    grid on
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % % Figure (5) Fluxes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure()
    subplot(2,1,1)
    errorbar(datenum(playa_filt.avg.t),mean(MOST_H,2,'omitnan'),...
        std(MOST_H,0,2,'omitnan'),'r-o')
    hold on
    errorbar(datenum(playa_filt.avg.t),mean(Shao_H,2,'omitnan'),...
        std(Shao_H,0,2,'omitnan'),'b-o')
    errorbar(datenum(playa_filt.avg.t),...
        mean(playa_filt.avg.H,2,'omitnan'),std(playa_filt.avg.H,0,2,'omitnan'),'k-o')
    legend('MOST','Shao','Tower')
    axis tight
    datetick('x','mmmdd HH:MM','keeplimits')
    ylabel('H (Wm$^{-2}$)')
    grid on
    
    
    subplot(2,1,2)
    errorbar(mean(MOST_H,1,'omitnan'),z,std(MOST_H,0,1,'omitnan'),'r-o','horizontal')
    hold on
    errorbar(mean(Shao_H,1,'omitnan'),z,std(Shao_H,0,1,'omitnan'),'b-o','horizontal')
    errorbar(mean(playa_filt.avg.H,1,'omitnan'),z,std(playa_filt.avg.H,0,1,'omitnan'),'k-o','horizontal')
    legend('MOST','Shao','Tower')
    ylabel('z (m)')
    xlabel('H (Wm$^{-2}$)')
    grid on
    axis tight
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
    if strcmp(day_of_interest,day_array(1))
        Shao_all = Shao_H;
        Most_all = MOST_H;
        wT_all = playa_filt.avg.H;
    else
        Shao_all = [Shao_all; Shao_H];
        Most_all = [Most_all; MOST_H];
        wT_all = [wT_all; playa_filt.avg.H];
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure (6) overall 1 to 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin = -10;xmax = 300;
ymin=xmin;ymax=xmax;
figure()
% errorbar(mean(wT_all,2,'omitnan'),mean(Shao_all,2,'omitnan'),std(Shao_all,0,2,'omitnan'),std(Shao_all,0,2,'omitnan'),std(wT_all,0,2,'omitnan'),std(wT_all,0,2,'omitnan'),'ro')
% hold on
% errorbar(mean(wT_all,2,'omitnan'),mean(Most_all,2,'omitnan'),std(Most_all,0,2,'omitnan'),std(Most_all,0,2,'omitnan'),std(wT_all,0,2,'omitnan'),std(wT_all,0,2,'omitnan'),'ko')
plot(mean(wT_all,2,'omitnan'),mean(Shao_all,2,'omitnan'),'ro')
hold on
plot(mean(wT_all,2,'omitnan'),mean(Most_all,2,'omitnan'),'ko')
% plot(mean(wT_all,2,'omitnan'),Shao_all(:,tower_height),'ro')
% hold on
% plot(mean(wT_all,2,'omitnan'),Most_all(:,tower_height),'ko')
plot(linspace(xmin, xmax, 100),linspace(xmin,xmax,100),'k-')
grid on
ylabel('Modeled $H$ (Wm$^{-2}$)')
xlabel('$\rho c_p \overline{w^\prime \theta^\prime}$ (Wm$^{-2}$)')
legend('Shao $\it{et}$ $\it{al.}$, n = 1.3','MOST')
axis([xmin xmax ymin ymax])
title('5 min')

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure (5) look at MOST Traditionally (Marc's paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_std = reshape(playa_filt.min_30.w, [length(playa_filt.min_30.w),]);

figure()
plot(std(playa_filt.min_30.w))








%% Look at Shao's dependance on stability
% % Compute the eddy diffusivity
% Kh_sg = -wTh_prime(:,5)./dTdz(5,:)';
% 
% 
% % Examine Shao's dependance on stability
% figure()
% subplot(1,2,1)
% plot(sqrt(playa_filt.min_30.tke(start_i:end_i-1,tower_height)).*z(tower_height),...
%     Kh_sg,'ko')
% hold on
% plot(linspace(.6,3, 100),linspace(.6,3,100),'k-')
% xlabel('$\sqrt{e} l$')
% ylabel('$K_h$')
% grid on
% axis([.6 3 .6 3])
% 
% subplot(1,2,2)
% plot(z(tower_height)./playa_filt.min_30.L(start_i:end_i-1,tower_height),...
%     Kh_sg./(sqrt(playa_filt.min_30.tke(start_i:end_i-1,tower_height)).*z(tower_height)),'ko')
% hold on
% plot(linspace(-4,0, 100),linspace(.4,.4,100),'k-')
% xlabel('$z/L$')
% ylabel('$\zeta = K_h/ \sqrt{e}l$')
% grid on
% axis tight


