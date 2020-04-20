clc
clear
close all

%%
directory_name = dir('H:\Montazeri\05 Unit_MRI_Data');
for k = 1:length(directory_name)
    if size(directory_name(k).name,2) == 21
        fname = directory_name(k).name;
        load(fname)
    end
end
clear directory_name fname k

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii = 1601;
ii = 801;

iddata_FRF_ACT_1L_SEN_1L = FRF_ACT_1L_SEN_1L(1:ii,3)'+FRF_ACT_1L_SEN_1L(1:ii,4)'*1i;
iddata_FRF_ACT_1L_SEN_2L = FRF_ACT_1L_SEN_2L(1:ii,3)'+FRF_ACT_1L_SEN_2L(1:ii,4)'*1i;
iddata_FRF_ACT_1L_SEN_3L = FRF_ACT_1L_SEN_3L(1:ii,3)'+FRF_ACT_1L_SEN_3L(1:ii,4)'*1i;
iddata_FRF_ACT_1L_SEN_1R = FRF_ACT_1L_SEN_1R(1:ii,3)'+FRF_ACT_1L_SEN_1R(1:ii,4)'*1i;
iddata_FRF_ACT_1L_SEN_2R = FRF_ACT_1L_SEN_2R(1:ii,3)'+FRF_ACT_1L_SEN_2R(1:ii,4)'*1i;
iddata_FRF_ACT_1L_SEN_3R = FRF_ACT_1L_SEN_3R(1:ii,3)'+FRF_ACT_1L_SEN_3R(1:ii,4)'*1i;

iddata_FRF_ACT_2L_SEN_1L = FRF_ACT_2L_SEN_1L(1:ii,3)'+FRF_ACT_2L_SEN_1L(1:ii,4)'*1i;
iddata_FRF_ACT_2L_SEN_2L = FRF_ACT_2L_SEN_2L(1:ii,3)'+FRF_ACT_2L_SEN_2L(1:ii,4)'*1i;
iddata_FRF_ACT_2L_SEN_3L = FRF_ACT_2L_SEN_3L(1:ii,3)'+FRF_ACT_2L_SEN_3L(1:ii,4)'*1i;
iddata_FRF_ACT_2L_SEN_1R = FRF_ACT_2L_SEN_1R(1:ii,3)'+FRF_ACT_2L_SEN_1R(1:ii,4)'*1i;
iddata_FRF_ACT_2L_SEN_2R = FRF_ACT_2L_SEN_2R(1:ii,3)'+FRF_ACT_2L_SEN_2R(1:ii,4)'*1i;
iddata_FRF_ACT_2L_SEN_3R = FRF_ACT_2L_SEN_3R(1:ii,3)'+FRF_ACT_2L_SEN_3R(1:ii,4)'*1i;

iddata_FRF_ACT_3L_SEN_1L = FRF_ACT_3L_SEN_1L(1:ii,3)'+FRF_ACT_3L_SEN_1L(1:ii,4)'*1i;
iddata_FRF_ACT_3L_SEN_2L = FRF_ACT_3L_SEN_2L(1:ii,3)'+FRF_ACT_3L_SEN_2L(1:ii,4)'*1i;
iddata_FRF_ACT_3L_SEN_3L = FRF_ACT_3L_SEN_3L(1:ii,3)'+FRF_ACT_3L_SEN_3L(1:ii,4)'*1i;
iddata_FRF_ACT_3L_SEN_1R = FRF_ACT_3L_SEN_1R(1:ii,3)'+FRF_ACT_3L_SEN_1R(1:ii,4)'*1i;
iddata_FRF_ACT_3L_SEN_2R = FRF_ACT_3L_SEN_2R(1:ii,3)'+FRF_ACT_3L_SEN_2R(1:ii,4)'*1i;
iddata_FRF_ACT_3L_SEN_3R = FRF_ACT_3L_SEN_3R(1:ii,3)'+FRF_ACT_3L_SEN_3R(1:ii,4)'*1i;

iddata_FRF_ACT_1R_SEN_1L = FRF_ACT_1R_SEN_1L(1:ii,3)'+FRF_ACT_1R_SEN_1L(1:ii,4)'*1i;
iddata_FRF_ACT_1R_SEN_2L = FRF_ACT_1R_SEN_2L(1:ii,3)'+FRF_ACT_1R_SEN_2L(1:ii,4)'*1i;
iddata_FRF_ACT_1R_SEN_3L = FRF_ACT_1R_SEN_3L(1:ii,3)'+FRF_ACT_1R_SEN_3L(1:ii,4)'*1i;
iddata_FRF_ACT_1R_SEN_1R = FRF_ACT_1R_SEN_1R(1:ii,3)'+FRF_ACT_1R_SEN_1R(1:ii,4)'*1i;
iddata_FRF_ACT_1R_SEN_2R = FRF_ACT_1R_SEN_2R(1:ii,3)'+FRF_ACT_1R_SEN_2R(1:ii,4)'*1i;
iddata_FRF_ACT_1R_SEN_3R = FRF_ACT_1R_SEN_3R(1:ii,3)'+FRF_ACT_1R_SEN_3R(1:ii,4)'*1i;

iddata_FRF_ACT_2R_SEN_1L = FRF_ACT_2R_SEN_1L(1:ii,3)'+FRF_ACT_2R_SEN_1L(1:ii,4)'*1i;
iddata_FRF_ACT_2R_SEN_2L = FRF_ACT_2R_SEN_2L(1:ii,3)'+FRF_ACT_2R_SEN_2L(1:ii,4)'*1i;
iddata_FRF_ACT_2R_SEN_3L = FRF_ACT_2R_SEN_3L(1:ii,3)'+FRF_ACT_2R_SEN_3L(1:ii,4)'*1i;
iddata_FRF_ACT_2R_SEN_1R = FRF_ACT_2R_SEN_1R(1:ii,3)'+FRF_ACT_2R_SEN_1R(1:ii,4)'*1i;
iddata_FRF_ACT_2R_SEN_2R = FRF_ACT_2R_SEN_2R(1:ii,3)'+FRF_ACT_2R_SEN_2R(1:ii,4)'*1i;
iddata_FRF_ACT_2R_SEN_3R = FRF_ACT_2R_SEN_3R(1:ii,3)'+FRF_ACT_2R_SEN_3R(1:ii,4)'*1i;

iddata_FRF_ACT_3R_SEN_1L = FRF_ACT_3R_SEN_1L(1:ii,3)'+FRF_ACT_3R_SEN_1L(1:ii,4)'*1i;
iddata_FRF_ACT_3R_SEN_2L = FRF_ACT_3R_SEN_2L(1:ii,3)'+FRF_ACT_3R_SEN_2L(1:ii,4)'*1i;
iddata_FRF_ACT_3R_SEN_3L = FRF_ACT_3R_SEN_3L(1:ii,3)'+FRF_ACT_3R_SEN_3L(1:ii,4)'*1i;
iddata_FRF_ACT_3R_SEN_1R = FRF_ACT_3R_SEN_1R(1:ii,3)'+FRF_ACT_3R_SEN_1R(1:ii,4)'*1i;
iddata_FRF_ACT_3R_SEN_2R = FRF_ACT_3R_SEN_2R(1:ii,3)'+FRF_ACT_3R_SEN_2R(1:ii,4)'*1i;
iddata_FRF_ACT_3R_SEN_3R = FRF_ACT_3R_SEN_3R(1:ii,3)'+FRF_ACT_3R_SEN_3R(1:ii,4)'*1i;

iddata_FRF_SHAKER_SEN_1L = FRF_SHAKER_SEN_1L(1:ii,3)'+FRF_SHAKER_SEN_1L(1:ii,4)'*1i;
iddata_FRF_SHAKER_SEN_2L = FRF_SHAKER_SEN_2L(1:ii,3)'+FRF_SHAKER_SEN_2L(1:ii,4)'*1i;
iddata_FRF_SHAKER_SEN_3L = FRF_SHAKER_SEN_3L(1:ii,3)'+FRF_SHAKER_SEN_3L(1:ii,4)'*1i;
iddata_FRF_SHAKER_SEN_1R = FRF_SHAKER_SEN_1R(1:ii,3)'+FRF_SHAKER_SEN_1R(1:ii,4)'*1i;
iddata_FRF_SHAKER_SEN_2R = FRF_SHAKER_SEN_2R(1:ii,3)'+FRF_SHAKER_SEN_2R(1:ii,4)'*1i;
iddata_FRF_SHAKER_SEN_3R = FRF_SHAKER_SEN_3R(1:ii,3)'+FRF_SHAKER_SEN_3R(1:ii,4)'*1i;

freq_vec = FRF_ACT_1L_SEN_1L(1:ii,2)*2*pi;
fmax = freq_vec(end);

iddata_ACT_1L = [iddata_FRF_ACT_1L_SEN_1L;iddata_FRF_ACT_1L_SEN_2L;iddata_FRF_ACT_1L_SEN_3L;iddata_FRF_ACT_1L_SEN_1R;iddata_FRF_ACT_1L_SEN_2R;iddata_FRF_ACT_1L_SEN_3R];
iddata_ACT_2L = [iddata_FRF_ACT_2L_SEN_1L;iddata_FRF_ACT_2L_SEN_2L;iddata_FRF_ACT_2L_SEN_3L;iddata_FRF_ACT_2L_SEN_1R;iddata_FRF_ACT_2L_SEN_2R;iddata_FRF_ACT_2L_SEN_3R];
iddata_ACT_3L = [iddata_FRF_ACT_3L_SEN_1L;iddata_FRF_ACT_3L_SEN_2L;iddata_FRF_ACT_3L_SEN_3L;iddata_FRF_ACT_3L_SEN_1R;iddata_FRF_ACT_3L_SEN_2R;iddata_FRF_ACT_3L_SEN_3R];
iddata_ACT_1R = [iddata_FRF_ACT_1R_SEN_1L;iddata_FRF_ACT_1R_SEN_2L;iddata_FRF_ACT_1R_SEN_3L;iddata_FRF_ACT_1R_SEN_1R;iddata_FRF_ACT_1R_SEN_2R;iddata_FRF_ACT_1R_SEN_3R];
iddata_ACT_2R = [iddata_FRF_ACT_2R_SEN_1L;iddata_FRF_ACT_2R_SEN_2L;iddata_FRF_ACT_2R_SEN_3L;iddata_FRF_ACT_2R_SEN_1R;iddata_FRF_ACT_2R_SEN_2R;iddata_FRF_ACT_2R_SEN_3R];
iddata_ACT_3R = [iddata_FRF_ACT_3R_SEN_1L;iddata_FRF_ACT_3R_SEN_2L;iddata_FRF_ACT_3R_SEN_3L;iddata_FRF_ACT_3R_SEN_1R;iddata_FRF_ACT_3R_SEN_2R;iddata_FRF_ACT_3R_SEN_3R];
iddata_SHAKER = [iddata_FRF_SHAKER_SEN_1L;iddata_FRF_SHAKER_SEN_2L;iddata_FRF_SHAKER_SEN_3L;iddata_FRF_SHAKER_SEN_1R;iddata_FRF_SHAKER_SEN_2R;iddata_FRF_SHAKER_SEN_3R];

%%
ny = 6;
nu = 7;
nf = ii;
N = nf;

Y = zeros(ny,nu,nf);

for i = 1:nf
    Y(:,1,i) = iddata_ACT_1L(:,i);
    Y(:,2,i) = iddata_ACT_2L(:,i);
    Y(:,3,i) = iddata_ACT_3L(:,i);
    Y(:,4,i) = iddata_ACT_1R(:,i);
    Y(:,5,i) = iddata_ACT_2R(:,i);
    Y(:,6,i) = iddata_ACT_3R(:,i);
    Y(:,7,i) = iddata_SHAKER(:,i);
end

w = real(freq_vec)';

n = 14;
m = nu;
p = ny;

T = 1/(100*2.56);

%%
z.y = Y; z.w = w;          % Specify measurements and frequencies they were obtained at
% z.T = T*0;       % Specify sample time

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how the estimation procedure runs
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ms.w = w; 
Ms.nx = n; 
Ms.op = 's'; 
Ms.T = T; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

identified_plant_ss = fsid(z,Ms);

%%
gt = identified_plant_ss;
A = gt.ss.A;
B = gt.ss.B(:,1:6);
H  = gt.ss.B(:,7);
C  = gt.ss.C;

[n,n_inp] = size(B);
[~,nw] = size(H);
[n_out,~] = size(C);

D = zeros(n_out,n_inp+nw);

gt_modal = canon(ss(A,[B H],C,D),'modal');
A = gt_modal.A;
B = gt_modal.B(:,1:6);
H  = gt_modal.B(:,7);
C  = gt_modal.C;

damp(gt_modal)
rank(obsv(A,C))
rank(ctrb(A,B))

% save identified An Bn Cn Hn D


for i = 1:nu
    for j = 1:ny
        figure
        semilogy(w,abs(squeeze(Y(j,i,:))))
        hold on
        semilogy(w,abs(squeeze(freqresp(gt_modal(j,i),w,'rad/s') )))
        xlim([0 fmax])
        xlabel('Frequency (rad/sec)')
        ylabel('FRF (V/V)')
    end
end