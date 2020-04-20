function error = q_n_finder(qn)
q = qn(1);
% n = qn(2);
r = qn(2);
q = round(q);
n=20;
% n = round(n);
r = round(r);

% %% Parameters
% m = 1;  % number of input
% p = 1;  % number of output
% 
% %% FFT parameters
% N = 6400; % 50,100,200,400,800,1600,3200,6400
% nFFT = 2.56*N;
% 
% % Baseband
% fmax = 800;  % frequency span (Hz)
% fs = 2.56*fmax; % Nyquist frequency (Hz)
% 
% dt = 1/fs;  % sampling time (sec)
% df = fmax/nFFT;  % sampling frequency (Hz)
% 
% %% Input signal generation with window and overlap
% load('freq_vec.mat')
% load('iddata_act1.mat')
% load('iddata_act2.mat')
% load('iddata_shaker.mat')
% 
% iddata_act1 = iddata_act1';
% iddata_act2 = iddata_act2';
% iddata_shaker = iddata_shaker';
% 
% [M,~] = size(freq_vec); % M-> number of frequency samples
% 
% G = zeros(p,m,M);
% if ((m == 1)&&(p == 1))
%     G(1,1,:) = iddata_act1(1,:);
% elseif ((m == 1)&&(p == 2))
%     G(:,1,:) = iddata_act1(:,:);
% elseif ((m == 2)&&(p == 1))
%     for i = 1:M
%         G(1,1,i) = iddata_act1(1,i);
%         G(1,2,i) = iddata_shaker(1,i);
%     end
% else
%     for i = 1:M
%     G(:,1,i) = iddata_act1(:,i);
%     G(:,2,i) = iddata_act2(:,i);
%     G(:,3,i) = iddata_shaker(:,i);
%     end
% end
% 
% %% Step 1: Discretization; frequency
% ww = (2*pi*freq_vec)';  % frequency vector from Modal analysis (rad/sec)
% ff = ww/(2*pi); % frequency vector from Modal analysis (Hz)
% 
% % transform continuous data to discrete using bilinear tansfrmation (Eq. (72))
% select_w = 1; % chon ferekans payin javabesh bade arc tan dar nazdike sefr javabesh khub nis.
% 
% % Please check
% % T = df*2*pi*0.02; % sort of like sampling period (Page 973)
% T = 5/max(ww);
% 
% wwd = 2*atan(T*ww(select_w:end)/2);
% wwd = wwd';
% Gz = G(:,:,select_w:end);
% Gzs = Gz;
% 
% % Eq. (11): Changing the discrete frequency range from [0 pi] to [0 2*pi]
% for k = 1:M-1
%     Gzs(:,:,M+k) = conj(Gz(:,:,M-k));
% end
% 
% %% Step 2: Construct h_hat in Eq. (12): 2M-point IDFT
% h_hat = zeros(p,m,2*M-1);
% for i = 1:2*M-1
%     h_hat_dum = zeros(p,m);
%     for k = 1:2*M-1
%         h_hat_dum = h_hat_dum + Gzs(:,:,k)*exp(1i*2*pi*i*k/(2*M));
%     end
%     h_hat(:,:,i) = h_hat_dum/(2*M);
%     i
% end
load matlab
%% Step 3 compute G_hat 

% Eq. (13) 
if q+r > 2*M
    error('q+r <= 2*M is not satisfied select smaller q and/or r (Eq. (13))')
end

H_hat = zeros(q*p,r*m);
for i = 1:q
    k = i;
    for j = 1:r
        H_hat((i-1)*p+1:i*p,(j-1)*m+1:j*m) = squeeze(h_hat(:,:,k));
        k = k+1;
    end
end

%% Step 4: SVD of H_hat
[U_hat,Sigma_hat,V_hat] = svd(H_hat);

%% Step 5: taking the n largest singular values
% According to equation (14) 
U_hat_s = U_hat(:,1:n);

%% Step 6: detemine A_hat & C_hat 
J1 = [eye((q-1)*p) zeros((q-1)*p,p)];
J2 = [zeros((q-1)*p,p) eye((q-1)*p)];
J3 = [eye(p) zeros(p,p*(q-1))];
A_hat = J1*U_hat_s;
A_hat = ((A_hat'*A_hat)^(-1))*A_hat';

A_hat = A_hat*J2*U_hat_s;
C_hat = J3*U_hat_s;

%% Step 7: detemine B_hat & D_hat  
% According to Eq. (20) solving batch LS problem 
S = [];
II = [];
for i = 1:(M-select_w+1)
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
for i = 1:(M-select_w+1)
    for j = 1:m
        BBB((k-1)*p+1:k*p,1) = squeeze(Gzs(:,j,i));
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

G_hat=ss(A_hat,B_hat,C_hat,D_hat,T);
sysc=d2c(G_hat,'tustin');
k = 1;
for i = 1:m
    for j = 1:p
        [magG,~]=freqresp(sysc(j,i),ff(select_w:end),'Hz');
        sss(:,k) = abs(squeeze(G(j,i,:)))-abs((squeeze(magG)));
        k = k+1;
%         figure;
%         semilogy(abs((squeeze(magG))))
%         hold on
%         semilogy(abs(squeeze(Gjw(j,i,select_w:end))))
    end
end
if (q < n+3) || (r < n+3)
    penalti = 1e20;
else
    penalti = 0;
end
error = norm(sss)+penalti;
