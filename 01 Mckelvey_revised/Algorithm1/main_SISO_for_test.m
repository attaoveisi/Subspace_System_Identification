clc
clear
% close all
load forsiso

%% SISO ID
ww=freq_vec*2*pi;
ff=ww/(2*pi);
j=sqrt(-1); 
H_2 = iddata_act1';
B_Complex=H_2;  %if acceleratoin is measured (FRF system az inja shuru mikonam man)
Gjw=B_Complex;   %if displacement is measured (chon acc andaze giri shode vali man lazem nadaram)
% =========================================================================
% transform continuous data to discrete using bilinear tansfrmation
% !!!!!!!!!!!!!!!! see help bilinear
select_w=1; % chon ferekans payin javabesh bade arc tan dar nazdike sefr javabesh khub nis.
T=5/max(ww);
wwd=2*atan(T*ww(select_w:end)/2);

Gzks=Gjw(select_w:end);
% Gzks(1,1)=500+500j;      % because Gjw(1,1) is inf         
% =========================================================================
m=1;                 % number of input
p=1;                 % number of output   
[M,MM]=size(Gzks);    % number of data
j=sqrt(-1);
n=16;                % order of system
q=n*5;        % q should be greater than system order (Eq 6)


% ===================== step 1: compute matrix G  =========================
% **************  according to equation (47) or (10) in conf paper*********
% G
for i=1:q
    for k=1:M
    Y_qp_mM(p*(i-1)+1:p*i,m*(k-1)+1:m*k)=exp(j*(i-1)*wwd(k,1))*Gzks((k-1)*p+1:k*p,:);
    end
end

%% ==================== step 1: compute matrix Wm  =========================
% *********************  according to equation (48) ***********************
% Wm
for i=1:q
    for k=1:M
    U_m(p*(i-1)+1:p*i,m*(k-1)+1:m*k)=exp(j*(i-1)*wwd(k,1))*eye(m,m);
    end
end
U_re=[real(U_m),imag(U_m)];
Y_re=[real(Y_qp_mM),imag(Y_qp_mM)];
% ================= step 2: compute QR factorizatin   =====================
% ******* according to equation (62) using lqdecomposition function *******
[Q,R]=qr([U_re' Y_re'],0);
R22=R(end-q+1:end,end-q+1:end);
% ==================== step 3: compute SVD  ===============================
% ******* according to equation (63) using lqdecomposition function *******
[U_hat,S_hat,V_hat] = svd(R22');

% =================== step 4: detemine system order  ======================
% **** according to equation (64) and estimate of observability matrix ****
U_hat_s=U_hat(:,1:n);

% ================= step 5: detemine A_hat & C_hat  =======================
% ***************** according to equations (65) & (66) ********************
O_hat_upperline=U_hat_s(p+1:end,:);
O_hat_underline=U_hat_s(1:end-p,:);
A_hat=inv(O_hat_underline'*O_hat_underline)*O_hat_underline'*O_hat_upperline;
C_hat=U_hat_s(1:p,:);
% ================= step 6: detemine B_hat & D_hat  =======================
% ********************* according to equation (67) ************************
teta(1:n+1,1)=0;       % the vector contains of B_hat and D_hat coefficient Az least square problem
for i=1:(M)
    C_hattimeinvA_hat(i,:)=C_hat*inv((exp(j*wwd(i,1)))*eye(size(A_hat))-A_hat);
end
zaribD(1:M,1)=1;
B_hat=teta(1:n,1);
D_hat=teta(n+1,1);
epsilon=Gzks(1:M)-D_hat*C_hattimeinvA_hat*B_hat; % estimation error
norm(epsilon)
jakoobiyan=[C_hattimeinvA_hat zaribD];
jakoobiyan_re=[real(jakoobiyan);imag(jakoobiyan)];
% teta=inv((jakoobiyan_re'*jakoobiyan_re))*((jakoobiyan_re')*[real(Gzks(1:M));imag(Gzks(1:M))]);
teta=inv((jakoobiyan_re'*jakoobiyan_re))*((jakoobiyan_re')*[real(Gzks(1:M));imag(Gzks(1:M))]);
B_hat=teta(1:n,1);
D_hat=teta(n+1,1);
epsilon=Gzks(1:M)-D_hat*C_hattimeinvA_hat*B_hat;
norm(epsilon)
% =============== step 7: estimated transfer function  ====================
% ********************* according to equation (68) ************************
G_hat=ss(A_hat,B_hat,C_hat,D_hat,T);
% figure
% bode(G_hat)

% % % % % % % % % % % % % % % % % % % % % %   
sysc=d2c(G_hat,'tustin');
figure
plot(ff(select_w:end),20*log10(abs(Gjw(select_w:end))));
hold on
[magG,phaseG,w]=bode(sysc,10.5*2*pi:pi:2*pi*800);
  for k=1:size(w)
      bb(k,1)=magG(:,:,k);
  end
  plot(w/2/pi,20*log10(bb),'r')
  
%   xlim([0 800])


