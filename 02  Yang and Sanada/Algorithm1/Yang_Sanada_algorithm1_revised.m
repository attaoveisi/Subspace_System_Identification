%% ===========================================================================================
% MIMO w-operator frequeny-domain SysID based on Algorithm 1 in Yang and Sanada 2000 
% in Elect. Eng. in Japan. The Script is written by ATTA OVEISI at Ruhr-Uni Bochum MAS in 
% Germany. Use of this code without permission is not allowed! Contact: atta.oveisi@rub.de
% ============================================================================================

%%
clc
clear
close all

%% Fast Fourrier Transformation parameters
% Baseband
fmax = 800;  % frequency span (Hz)
fs = 2.56*fmax; % Nyquist frequency (Hz)

dt = 1/fs;  % sampling time (sec)

m = 1;                 % number of input
l = 1;                 % number of output   

% Define ii & alfa & n
n = 12;
ii = 36;
alfa = 800*2*pi; % 5
beta = 47;
if ii <= n
    error('ii should be greater than n');
end

if beta < n/m
    error('beta must be greater than n/m')
end

%% Input signal generation with window and overlap
load('iddata.mat')
freq_vec = act1_1(:,2);

ww = freq_vec*2*pi;
T = 5/max(ww);
select_w = 1;
wk = 2*atan(T*ww(select_w:end)/2);

iddata_act1 = (act1_1(:,3)+1i*act1_1(:,4))';
iddata_act2 = (act2_1(:,3)+1i*act2_1(:,4))';
iddata_shaker = (shaker_1(:,3)+1i*shaker_1(:,4))';

% Up to freq(index) Hz 
index = 6401;
freq_vec = freq_vec(1:index,:);
iddata_act1 = iddata_act1(:,1:index);
iddata_act2 = iddata_act2(:,1:index);
iddata_shaker = iddata_shaker(:,1:index);

[M,~] = size(freq_vec); % M-> number of frequency samples

Gm = zeros(l,m,M);
if ((m == 1)&&(l == 1))
    Gm(1,1,:) = iddata_act1(1,:);
elseif ((m == 1)&&(l == 2))
    Gm(:,1,:) = iddata_act1(:,:);
elseif ((m == 2)&&(l == 1))
    for i = 1:M
        Gm(1,1,i) = iddata_act1(1,i);
        Gm(1,2,i) = iddata_shaker(1,i);
    end
else
    for i = 1:M
    Gm(:,1,i) = iddata_act1(:,i);
    Gm(:,2,i) = iddata_act2(:,i);
    Gm(:,3,i) = iddata_shaker(:,i);
    end
end

%% Step (1)
%Construct G_i,M^f (Eq. (19))
Gm_f = zeros(size(Gm));
Im = ones(m,m);
Im_f = zeros(m,M*m);

omega = (freq_vec*2*pi)';

% ww = freq_vec*2*pi;
% T = 5/max(ww);
% select_w = 1;
% omega = 2*atan(T*ww(select_w:end)/2)';

omega_alfa = @(omegai,alfai) (1i*omegai-alfai)/(1i*omegai+alfai); 

for i = 1:m
    for j = 1:l
        for k = 1:M
            Gm_f(j,i,k) = (sqrt(2)*alfa/(1i*omega(1,k)+alfa))*Gm(j,i,k); % Eq. (17)
            if i == 1 && j == 1
                Im_f(:,(k-1)*m+1:k*m) = (sqrt(2)*alfa/(1i*omega(1,k)+alfa))*Im;
            end
        end
    end
end

I_iM_f = zeros(ii*m,M*m);
G_iM_f = zeros(ii*l,M*m);
for i = 1:ii
    for j = 1:M
        G_iM_f((i-1)*l+1:i*l,(j-1)*m+1:j*m) = ((omega_alfa(omega(1,j),alfa))^(i-1))*squeeze(Gm_f(:,:,j));
        I_iM_f((i-1)*m+1:i*m,(j-1)*m+1:j*m) = ((omega_alfa(omega(1,j),alfa))^(i-1))*Im_f(:,(j-1)*m+1:j*m);
    end
end
G_iM_f = [real(G_iM_f) imag(G_iM_f)]; % Eq. (24)
I_iM_f = [real(I_iM_f) imag(I_iM_f)]; % Eq. (24)

%% Step (2): LQ factorization 
% Refer to Eq. (27)

% !!!!!!!! Wrong: This must be LQ factorization
[Q_dum,R_dum] = qr([I_iM_f;G_iM_f]',0); % !!!!!!!! Wrong
R_dum = R_dum';
R_11 = R_dum(1:m*ii,1:m*ii);
R_12 = R_dum(1:m*ii,m*ii+1:end);
R_21 = R_dum(m*ii+1:end,1:m*ii);
R_22 = R_dum(m*ii+1:end,m*ii+1:end);

%% Step (3): SVD of R_22
% Eq. (32)
[U_n,S_n,V_n] = svd(R_22);

V_n_t = V_n';

Un = U_n(:,1:n);
Un_p = U_n(:,n+1:end);

Sn = S_n(1:n,1:n);
S2 = S_n(n+1:end,n+1:end);

Vn_t = V_n_t(1:n,:);
Vn_t_p = V_n_t(n+1:end,:);
Vn = Vn_t';
Vn_p = Vn_t_p';

%% Step (4): Construct Un(1) & Un(2) and calculated A_wT & C_wT
% Eqs. (34) and (35)
Un1 = Un(1:(ii-1)*l,:);
Un2 = Un(end-(ii-1)*l+1:end,:);

X_plus = @(X) ((X'*X)^(-1))*X'; % Moore-Penrose pseudoinverse of full-column rand matrix X

A_wT = (X_plus(Un1))*Un2;
C_wT = Un(1:l,:);

%% Step (5): Compute B_wT and D_wT
Li = zeros(M*l,m);
for i = 1:M
    % Eq. (37)
    Li((i-1)*l+1:i*l,:) = squeeze(Gm(:,:,i));
end

% For strictly proper system: D_c = 0:
Mi2 = zeros(l*M,n);
for i = 1:M
    % Eq. (41)
    Mi2((i-1)*l+1:i*l,:) = C_wT*(inv((omega_alfa(omega(1,i),alfa))*eye(n)-A_wT)-inv(eye(n)-A_wT));
end

% (a) Construct the matrices inside Eq. (39)
LHS = [real(Li);imag(Li)];
RHS = [real(Mi2);imag(Mi2)];

B_wT = (X_plus(RHS))*LHS;
D_wT = -(C_wT/(eye(n)-A_wT))*B_wT;

%% Step (6): Continuous system
G_hat = ss(A_wT,B_wT,C_wT,D_wT,T);
% sysc = d2c(G_hat,'tustin');
% figure(2)
% bode(G_hat)
% hold on
% bode(sysc)
A_c = (alfa*((eye(n)-A_wT)^(-1))*(eye(n)+A_wT));
B_c = (sqrt(2*alfa)*((eye(n)-A_wT)^(-1))*B_wT);
C_c = (sqrt(2*alfa)*C_wT*((eye(n)-A_wT)^(-1)));
D_c = (D_wT+C_wT*((eye(n)-A_wT)^(-1))*B_wT);
sysc = ss(A_c,B_c,C_c,D_c);

figure(3)
k = 0;
for i = 1:m
    for j = 1:l
        k = k+1;
        subplot(m,l,k)
        plot(freq_vec,20*log10(abs(squeeze(Gm(j,i,1:M)))));
        hold on
        [magG,phaseG,w] = bode(sysc(j,i),freq_vec*2*pi);
        for jj = 1:size(w)
            bb(jj,1) = squeeze(magG(:,:,jj));
        end
        plot(w/2/pi,20*log10(bb),'r')
        xlim([0 fmax])
    end
end