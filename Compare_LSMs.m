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
load('.\Materhorn_data\playaSpring30minLinDetUTESpac3.mat');
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
        error_shao(i,z) = ((H_tower(i,z) - mean(result_shao(z).H(start_index:end_index)))/H_tower(i,z))*100;
        error_MOST(i,z) = ((H_tower(i,z) - mean(result_MOST(z).H(start_index:end_index)))/H_tower(i,z))*100;
        start_index = end_index+1;
        end_index = end_index+30;
    end
end



%%
figure()
errorbar(mean(error_shao(:,1),2),std(error_shao'))
hold on 
errorbar(mean(error_MOST(:,1),2),std(error_MOST'))



%%


figure()
boxplot(error_shao')
hond on 
boxplot(error_MOST')
