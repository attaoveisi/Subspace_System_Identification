clc
clear
close all

%% Fast Fourrier Transformation parameters
N = 6400; % 50,100,200,400,800,1600,3200,6400
nFFT = 2.56*N;

% Baseband
fmax = 800;  % frequency span (Hz)
fs = 2.56*fmax; % Nyquist frequency (Hz)

dt = 1/fs;  % sampling time (sec)
df = fmax/nFFT;  % sampling frequency (Hz)

m=1;                 % number of input
p=2;                 % number of output   

%% Input signal generation with window and overlap
load('freq_vec.mat')
load('iddata_act1.mat')
load('iddata_act2.mat')
load('iddata_shaker.mat')

iddata_act1 = iddata_act1';
iddata_act2 = iddata_act2';
iddata_shaker = iddata_shaker';

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


%% MIMO ID
ww = (2*pi*freq_vec)'; 
ff = ww/(2*pi);
Gjw = G;

%% =========================================================================
% transform continuous data to discrete using bilinear tansfrmation
% !!!!!!!!!!!!!!!! see help bilinear
select_w = 50; % chon ferekans payin javabesh bade arc tan dar nazdike sefr javabesh khub nis.
T = 5/max(ww);
wwd = 2*atan(T*ww(select_w:end)/2);
wwd = wwd';
Gzks = Gjw(:,:,select_w:end);

%% =========================================================================
n = 16;                % order of system
q = 18;        % q should be greater than system order (Eq 6)

%% ===================== step 1: compute matrix G  =========================
% **************  according to equation (47) or (10) in conf paper*********
for i = 1:q
    for k = 1:M-select_w+1
    Y_qp_mM(p*(i-1)+1:p*i,m*(k-1)+1:m*k) = exp(1i*(i-1)*wwd(k,1))*squeeze(Gzks(:,:,k));
    end
end

%% ==================== step 1: compute matrix Wm  =========================
% *********************  according to equation (48) ***********************
for i=1:q
    for k=1:M-select_w+1
    U_m(m*(i-1)+1:m*i,m*(k-1)+1:m*k) = exp(1i*(i-1)*wwd(k,1))*eye(m,m);
    end
end
U_re = [real(U_m),imag(U_m)];
Y_re = [real(Y_qp_mM),imag(Y_qp_mM)];

%% ================= step 2: compute QR factorizatin   =====================
% ******* according to equation (62) using lqdecomposition function *******
[Q,R] = qr([U_re' Y_re'],0);
R22 = R(end-q+1:end,end-q+1:end);

%% ==================== step 3: compute SVD  ===============================
% ******* according to equation (63) using lqdecomposition function *******
[U_hat,S_hat,V_hat] = svd(R22');

%% =================== step 4: detemine system order  ======================
% **** according to equation (64) and estimate of observability matrix ****
U_hat_s = U_hat(:,1:n);

%% ================= step 5: detemine A_hat & C_hat  =======================
% ***************** according to equations (65) & (66) ********************
O_hat_upperline = U_hat_s(p+1:end,:);
O_hat_underline = U_hat_s(1:end-p,:);
A_hat = inv(O_hat_underline'*O_hat_underline)*O_hat_underline'*O_hat_upperline;
C_hat = U_hat_s(1:p,:);

%% ================= step 6: detemine B_hat & D_hat  =======================
% ********************* according to equation (67) ************************
S = [];
II = [];
for i = 1:M-select_w+1
    S = [S ; C_hat*inv((exp(1i*wwd(i,1)))*eye(size(A_hat))-A_hat)];
    II = [II ; eye(p)];
end
AA = [S II];
AAA = [];
for i = 1:m
    AAA = blkdiag(AAA,AA);
end

AAAA = [real(AAA);imag(AAA)];

k = 1;
for i = 1:M-select_w+1
    for j = 1:m
        BBB((k-1)*p+1:k*p,1) = squeeze(Gzks(:,j,i));
        k = k+1;
    end
end

BBBB = [real(BBB);imag(BBB)];

teta = ((AAAA'*AAAA)^-1)*AAAA'*BBBB;

for i = 1:m
    teta_dum = teta((i-1)*(p+n)+1:i*(p+n),1);
    B_hat(:,i) = teta_dum(1:n,:);
    D_hat(:,i) = teta_dum(n+1:end,:);
end

for i = 1:M-select_w+1
    epsilon(:,:,i) = squeeze(Gzks(:,:,i))-D_hat-(C_hat/((exp(1i*wwd(i,1)))*eye(size(A_hat))-A_hat))*B_hat; % estimation error
    normed_eps(1,i) = norm(squeeze(epsilon(:,:,i)));
end
figure(1);
plot(ff(select_w:end),normed_eps);

%% =============== step 7: estimated transfer function  ====================
% ********************* according to equation (68) ************************
G_hat=ss(A_hat,B_hat,C_hat,D_hat,T);
% figure
% bode(G_hat)

% % % % % % % % % % % % % % % % % % % % % %   
sysc=d2c(G_hat,'tustin');
figure(2)
k = 0;
for i = 1:m
    for j = 1:p
        k = k+1;
        subplot(m,p,k)
        plot(ff(select_w:end),20*log10(abs(squeeze(Gjw(j,i,select_w:end)))));
        hold on
        [magG,phaseG,w]=bode(sysc(j,i),ff(select_w:end)*2*pi);
        for l = 1:size(w)
            bb(l,1) = squeeze(magG(:,:,l));
        end
        plot(w/2/pi,20*log10(bb),'r')
        xlim([0 800])
    end
end