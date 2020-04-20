%%
clc
clear
close all

%% Fast Fourrier Transformation parameters
% Baseband
fmax = 800;  % frequency span (Hz)
fs = 2.56*fmax; % Nyquist frequency (Hz)

dt = 1/fs;  % sampling time (sec)

m = 2;                 % number of input
p = 1;                 % number of output   

%% Input signal generation with window and overlap
load('iddata.mat')
freq_vec = act1_1(:,2);
ww = freq_vec*2*pi;

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

G = zeros(p,m,M);
if ((m == 1)&&(p == 1))
    G(1,1,:) = iddata_act1(1,:);
elseif ((m == 1)&&(p == 2))
    G(:,1,:) = iddata_act1(:,:);
elseif ((m == 2)&&(p == 1))
    for i = 1:M
        G(1,1,i) = iddata_act1(1,i);
        G(1,2,i) = iddata_shaker(1,i);
    end
else
    for i = 1:M
    G(:,1,i) = iddata_act1(:,i);
    G(:,2,i) = iddata_act2(:,i);
    G(:,3,i) = iddata_shaker(:,i);
    end
end


%% Uniformly-spaced frequencies
T = 5/max(ww);
select_w = 1;
wk = 2*atan(T*ww(select_w:end)/2);

% Assign system parameters
n = 12;                % order of system
q = 5*n;        % q should be greater than system order (Eq 6)

%% ===================== step 1: compute matrix G  =========================
% **************  according to equation (47) or (10) in conf. paper*********
G_qp_mM = zeros(q*p,m*M);
for i = 1:q
    for k = 1:M
%         disp([i k])
        G_qp_mM(p*(i-1)+1:p*i,m*(k-1)+1:m*k) = exp(1i*(i-1)*wk(k,1))*squeeze(G(:,:,k));
    end
end
% G_qp_mM = G_qp_mM./sqrt(M);

%% ==================== step 2: compute matrix Wm  =========================
% *********************  according to equation (48) ***********************
W_m = zeros(q*m,m*M);
for i = 1:q
    for k = 1:M
        W_m(m*(i-1)+1:m*i,m*(k-1)+1:m*k) = exp(1i*(i-1)*wk(k,1))*eye(m,m);
    end
end
% W_m = W_m./sqrt(M);

%% ==================== step 3: Convert to real-valued equation  =========================
% *********************  according to equation (51) ***********************
Gi = [real(G_qp_mM) imag(G_qp_mM)];
Wi = [real(W_m) imag(W_m)];

%% ================= step 4: compute QR factorizatin   =====================
% ******* according to equation (62) using lq decomposition function *******
[Q_dum,R_dum] = qr([Wi' Gi'],0);
R_11 = R_dum(1:end-q,1:end-q);
R_12 = R_dum(1:end-q,end-q+1:end);
R_21 = R_dum(end-q+1:end,1:end-q);
R_22 = R_dum(end-q+1:end,end-q+1:end);

%% ==================== step 5: compute SVD  ===============================
% ******* according to equation (63) using lqdecomposition function *******
[U_hat,Sigma_hat,V_hat] = svd(R_22');

%% =================== step 6: detemine system order  ======================
% **** according to equation (64) and estimate of observability matrix ****
V_hat_t = V_hat';

U_hat_s = U_hat(:,1:n);
U_hat_o = U_hat(:,n+1:end);

Sigma_hat_s = Sigma_hat(1:n,1:n);
Sigma_hat_n = Sigma_hat(n+1:end,n+1:end);

V_hat_t_s = V_hat_t(1:n,:);
V_hat_t_o = V_hat_t(n+1:end,:);
V_hat_s = V_hat_t_s';
V_hat_o = V_hat_t_o';

%% ================= step 7: detemine A_hat & C_hat  =======================
% Eqs. (17)-(19)
J1 = [eye((q-1)*p) zeros((q-1)*p,p)];
J2 = [zeros((q-1)*p,p) eye((q-1)*p)];
J3 = [eye(p) zeros(p,(q-1)*p)];

% Eqs (15) and (16): Finding A_hat and C_hat
X_plus = @(X) ((X'*X)^(-1))*X'; % Moore-Penrose pseudoinverse of full-column rand matrix X

A_hat = (X_plus(J1*U_hat_s))*J2*U_hat_s;
C_hat = J3*U_hat_s;

%% ================= step 8: detemine B_hat & D_hat  =======================
% ********************* according to equation (67) ************************
S = [];
II = [];
for i = 1:M
    S = [S ; C_hat*inv((exp(1i*wk(i,1)))*eye(size(A_hat))-A_hat)];
    II = [II ; eye(p)];
end
AA = [S II];
AAA = [];

for i = 1:m
    AAA = blkdiag(AAA,AA);
end

AAAA = [real(AAA);imag(AAA)];

k = 1;
for i = 1:m
    for j = 1:M
        BBB((k-1)*p+1:k*p,1) = squeeze(G(:,i,j));
        k = k+1;
    end
end

BBBB = [real(BBB);imag(BBB)];

teta = ((AAAA'*AAAA)^-1)*AAAA'*BBBB;

B_hat = zeros(n,m);
D_hat = zeros(p,m);
for i = 1:m
    teta_dum = teta((i-1)*(p+n)+1:i*(p+n),1);
    B_hat(:,i) = teta_dum(1:n,:);
    D_hat(:,i) = teta_dum(n+1:end,:);
end

%% =============== step 7: estimated transfer function  ====================
% ********************* according to equation (68) ************************
G_hat = ss(A_hat,B_hat,C_hat,D_hat,T);
sysc = d2c(G_hat,'tustin');
figure(2)
bode(G_hat)
hold on
bode(sysc)

figure(3)
k = 0;
for i = 1:m
    for j = 1:p
        k = k+1;
        subplot(m,p,k)
        plot(freq_vec,20*log10(abs(squeeze(G(j,i,:)))));
        hold on
        [magG,phaseG,w] = bode(sysc(j,i),freq_vec*2*pi);
        for l = 1:size(w)
            bb(l,1) = squeeze(magG(:,:,l));
        end
        plot(w/2/pi,20*log10(bb),'r')
        xlim([0 fmax])
    end
end