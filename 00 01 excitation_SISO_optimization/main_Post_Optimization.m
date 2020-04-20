clc
clear
close all

%% 
load postotimization
load identified

%%
n = 4;
dt = 1e-4; 
total_time = 30;
a = zeros(1,14);
opts=simset('Solver','ode8','FixedStep',dt,'SrcWorkspace','current','DstWorkspace','current');
sim('Sim_Func_g',[0 total_time],opts);
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
plot(y(1,:),y(2,:))
hold on
plot(output_model.time,output_model.signals.values)
xlim([0 total_time])
legend('measurement','Frequency-domain SysID+Parameter identification')
%%
load postotimization

a = x';
% a = [556755.252099794,-52827305.5478599,0.211818369465259,46.3760343805078,24.9233402579705,-5.05951560717271,-2604.74248077992,34.8657393411266,-14.4323548096724,-8765223.40197508,-14755662.5900491,-2883.73134915326,229.095081260287,-11955.4677328915];
opts=simset('Solver','ode8','FixedStep',dt,'SrcWorkspace','current','DstWorkspace','current');
sim('Sim_Func_g',[0 total_time],opts);
subplot(2,1,2)
plot(y(1,:),y(2,:))
hold on
plot(output_model.time,output_model.signals.values)
xlim([0 total_time])
legend('measurement','...+Geometrical nonlinearity')