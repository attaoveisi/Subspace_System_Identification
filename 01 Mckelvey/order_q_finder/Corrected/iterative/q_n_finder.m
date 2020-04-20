function [select_w,ff,G,sysc] = q_n_finder(q,r)

load matlab
select_w=1;

%% Step 3 compute G_hat 
% Parameters
n = 20; % identified system order


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

%% Final Step: estimated transfer function 
G_hat=ss(A_hat,B_hat,C_hat,D_hat,T);
sysc=d2c(G_hat,'tustin');
