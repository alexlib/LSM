%Compare heights and error between shao and MOST models

clear; close all;clc;

driving_data_set = 'MATERHORN';
HF_option_shao = 'Shao';
HF_option_MOST = 'MOST';
plots = 'off';
atmrespondTF = 'on';
ABLTF= 'on';
t_day = 1;
tower_height = 1;
for i = 1:6
    [result_shao(i)] = LSM(driving_data_set,HF_option_shao,plots,atmrespondTF,ABLTF,t_day,i);
end

for i = 1:6
    [result_MOST(i)] = LSM(driving_data_set,HF_option_MOST,plots,atmrespondTF,ABLTF,t_day,i);
end


%Compute %error from obs data
load('Materhorn_data/playaSpring30minLinDetUTESpac3.mat');
%%  %load tower data: 30 min data
H_index = [15 27 75 39 51 63];
%15 for 25.5 m, 27 for 19.4 m, 63 for 0.61m, 75 for 10.4 m,39 for 5 m,
%51 for 2 m
data_start = 1072;%25 May 0000 UTC
data_end = data_start + (2*24);%26 May 0000 UTC
H_tower = playaSpring.H(data_start:data_end,2).*...
    playaSpring.H(data_start:data_end,3).*playaSpring.H(data_start:data_end,H_index);
for z = 1:6
    start_index = 1;
    end_index = 30;
    for i = 1:length(H_tower)-1
        error_shao(i,z) = (( mean(result_shao(z).H(start_index:end_index))-H_tower(i,z))/H_tower(i,z))*100;
        error_MOST(i,z) = (( mean(result_MOST(z).H(start_index:end_index))-H_tower(i,z))/H_tower(i,z))*100;
        start_index = end_index+1;
        end_index = end_index+30;
    end
end



%% Make plots
close all; 
for z= 1:6
    
heights = [25 19 10 5 2 0.61];
cnt = length(result_shao(1).T(1:end));
% formatting
addpath('/Users/travismorrison/Documents/Code/functions_library')
ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);

%plot Simulation Results Shao
figure(1)
subplot(2,2,1)
%title('Simulated Temperature')
plot(linspace(0,t_day,cnt),result_shao(z).T,'-k')
hold on
plot(linspace(0,t_day,cnt),result_shao(z).Ta,'--k')
legend('$T_s$','$T_a$')
ylabel('T [K]')


subplot(2,2,2)
%title('Simulated SHF and LHF')
plot(linspace(0,t_day,cnt),result_shao(z).H,'-r')
hold on
plot(linspace(0,t_day,cnt),result_shao(z).LH,'-b')
legend('H','LH')

ylabel('E [Wm$^{-2}$]')

subplot(2,2,3)
%title('Simulated ABL')
plot(linspace(0,t_day,cnt),result_shao(z).ABL_h,'--k')
hold on
plot(linspace(0,t_day,cnt),linspace(2500,2500,cnt),'-k')
legend('Simulated','Max Obs')
xlabel('time [days]')
ylabel('h [m]')
axis([0 t_day 0 3000])

subplot(2,2,4)
%title('Simulated CV budget (positive for into CV)')
plot(linspace(0,t_day,cnt),result_shao(z).Rnet,'-k')
hold on
plot(linspace(0,t_day,cnt),-result_shao(z).G,'-g')
plot(linspace(0,t_day,cnt),-result_shao(z).H,'-r')
plot(linspace(0,t_day,cnt),-result_shao(z).LH,'-b')
legend('$R_{net}$','G','H','$H_L$')
xlabel('time [days]')
ylabel('E [Wm$^{-2}$]')
axis tight

[ax]=subtitle(['Tower Height: ',num2str(heights(z)),' m, Model: Shao']);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);
saveas(figure(1),[pwd '/plots_results/',num2str(heights(z)),'m/Shao_model_output.fig']);
saveas(figure(1),[pwd '/plots_results/',num2str(heights(z)),'m/Shao_model_output.png']);
%
%plot Simulation Results Shao
figure(2)
subplot(2,2,1)
%title('Simulated Temperature')
plot(linspace(0,t_day,cnt),result_MOST(z).T,'-k')
hold on
plot(linspace(0,t_day,cnt),result_MOST(z).Ta,'--k')
legend('$T_s$','$T_a$')
ylabel('T [K]')


subplot(2,2,2)
%title('Simulated SHF and LHF')
plot(linspace(0,t_day,cnt),result_MOST(z).H,'-r')
hold on
plot(linspace(0,t_day,cnt),result_MOST(z).LH,'-b')
legend('H','LH')

ylabel('E [Wm$^{-2}$]')

subplot(2,2,3)
%title('Simulated ABL')
plot(linspace(0,t_day,cnt),result_MOST(z).ABL_h,'--k')
hold on
plot(linspace(0,t_day,cnt),linspace(2500,2500,cnt),'-k')
legend('Simulated','Max Obs')
xlabel('time [days]')
ylabel('h [m]')
axis([0 t_day 0 3000])

subplot(2,2,4)
%title('Simulated CV budget (positive for into CV)')
plot(linspace(0,t_day,cnt),result_MOST(z).Rnet,'-k')
hold on
plot(linspace(0,t_day,cnt),-result_MOST(z).G,'-g')
plot(linspace(0,t_day,cnt),-result_MOST(z).H,'-r')
plot(linspace(0,t_day,cnt),-result_MOST(z).LH,'-b')
legend('$R_{net}$','G','H','$H_L$')
xlabel('time [days]')
ylabel('E [Wm$^{-2}$]')
axis tight

[ax]=subtitle(['Tower Height: ',num2str(heights(z)),' m, Model: MOST']);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);
saveas(figure(1),[pwd '/plots_results/',num2str(heights(z)),'m/MOST_model_output.fig']);
saveas(figure(1),[pwd '/plots_results/',num2str(heights(z)),'m/MOST_model_output.png']);

% plot %error for SHF
figure(3)
semilogy(linspace(0,24,48), abs(error_shao(:,z)),'k*-')
hold on
semilogy(linspace(0,24,48),abs(error_MOST(:,z)),'bo-')
legend('shao','MOST')
axis tight
ylabel('$\%$ Error H')
xlabel('time [hour]')
title([num2str(heights(z)),' m '])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);
saveas(figure(3),[pwd '/plots_results/',num2str(heights(z)),'m/H_error.fig']);
saveas(figure(3),[pwd '/plots_results/',num2str(heights(z)),'m/H_error.png']);

%
%compare Shao against experimental
H_index = [15 27 75 39 51 63];
T_index =[ 2 4 6 8 10 12];

figure(4)
subplot(2,2,1)
plot(linspace(0,24,1440),result_shao(z).Ta,'-k')
hold on
plot(linspace(0,24.5,49),playaSpring.Playa_1HZ(data_start:data_end,T_index(z))+273.15,'--k')
legend('$T_a$ modeled','$T_a$ obs')
ylabel('T [K]')
axis tight

subplot(2,2,2)
plot(linspace(0,24,cnt/t_day), result_shao(z).H((cnt-cnt/t_day+1):cnt),'-k')
hold on
plot(linspace(0,24.5,49),playaSpring.H(data_start:data_end,2).*...
    playaSpring.H(data_start:data_end,3).*playaSpring.H(data_start:data_end,H_index(z)),'--k')
%15 for 25.5 m, 27 for 19.4 m, 63 for 0.61m, 75 for 10.4 m,39 for 5 m,
%51 for 2 m
legend('modeled','obs ')
axis tight
ylabel('H [Wm$^{-2}$]')
xlabel('time [hrs]')

subplot(2,2,3)
plot(linspace(0,24,cnt/t_day), result_shao(z).LH,'-k')
hold on
plot(linspace(0,24.5,49),playaSpring.LHflux(data_start:data_end,6),'--k')
legend('modeled','obs (10.4 m)')
axis tight
ylabel('$H_L$ [Wm$^{-2}$]')
xlabel('time [hrs]')
axis tight
[ax]=subtitle(['Tower Height: ',num2str(heights(z)),' m, Model: Shao']);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);

saveas(figure(4),[pwd '/plots_results/',num2str(heights(z)),'m/Shao_term_comparison.fig']);
saveas(figure(4),[pwd '/plots_results/',num2str(heights(z)),'m/Shao_term_comparison.png']);

%compare MOST against obs
figure(5)
subplot(2,2,1)
plot(linspace(0,24,1440),result_MOST(z).Ta,'-k')
hold on
plot(linspace(0,24.5,49),playaSpring.Playa_1HZ(data_start:data_end,T_index(z))+273.15,'--k')
legend('$T_a$ modeled','$T_a$ obs')
ylabel('T [K]')
axis tight

subplot(2,2,2)
plot(linspace(0,24,cnt/t_day), result_MOST(z).H((cnt-cnt/t_day+1):cnt),'-k')
hold on
plot(linspace(0,24.5,49),playaSpring.H(data_start:data_end,2).*...
    playaSpring.H(data_start:data_end,3).*playaSpring.H(data_start:data_end,H_index(z)),'--k')
%15 for 25.5 m, 27 for 19.4 m, 63 for 0.61m, 75 for 10.4 m,39 for 5 m,
%51 for 2 m
legend('modeled','obs ')
axis tight
ylabel('H [Wm$^{-2}$]')
xlabel('time [hrs]')

subplot(2,2,3)
plot(linspace(0,24,cnt/t_day), result_MOST(z).LH,'-k')
hold on
plot(linspace(0,24.5,49),playaSpring.LHflux(data_start:data_end,6),'--k')
legend('modeled','obs (10.4 m)')
axis tight
ylabel('$H_L$ [Wm$^{-2}$]')
xlabel('time [hrs]')
axis tight
[ax]=subtitle(['Tower Height: ',num2str(heights(z)),' m, Model: MOST']);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);

saveas(figure(5),[pwd '/plots_results/',num2str(heights(z)),'m/MOST_term_comparison.fig']);
saveas(figure(5),[pwd '/plots_results/',num2str(heights(z)),'m/MOST_term_comparison.png']);
close all
end
%% Plot mean error

hr_start = 7*2;
hr_end = 16*2;
figure()
handaxes1 = axes('Position', [0.12 0.12 0.8 0.8]);
plot(mean(abs(error_shao(hr_start:hr_end,:)),1),heights,'k-*')
hold on
plot(mean(abs(error_MOST(hr_start:hr_end,:)),1),heights,'k--o')
legend('Shao','MOST')
xlabel('$\%$ Error during convective period')
ylabel('Height (m)')
grid on
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .7, 0.96]);
% Place second set of axes on same plot
handaxes2 = axes('Position', [0.7 0.3 0.2 0.3]);

plot(linspace(0,24,1440),result_shao(1).SWdn,  'k')
hold on
X = [linspace(hr_start/2,hr_end/2,(hr_end-hr_start)*30),fliplr(linspace(hr_start/2,hr_end/2,(hr_end-hr_start)*30))];
Y = [linspace(0,0,(hr_end-hr_start)*30), fliplr(result_shao(1).SWdn(hr_start*30:hr_end*30-1))];
fill(X,Y,[.7 .7 .7])
%set(handaxes2, 'Box', 'off')
xlabel('t (hrs)')
ylabel('SW$\downarrow$ (Wm$^{-2}$)')
axis tight

%%
close all









