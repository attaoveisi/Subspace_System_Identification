function ID_error = opt_func(inp)
% Fast Fourrier Transformation parameters
% Baseband
q = ceil(inp(1,1));
r = ceil(inp(1,2));

m = 1;                 % number of input
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

% Uniformly-spaced frequencies
% wk = ((0:M)*pi/M)';
T = 5/max(ww);
select_w = 1;
wk = 2*atan(T*ww(select_w:end)/2);

% Eq. (11): Extend the transfer function to the full unit circle.
for k = 1:M-1
    
    for i = 1:p
        for j = 1:m
            G(i,j,M+k) = conj(G(i,j,M-k)); % if MIMO maybe a transpose is required
        end
    end
end

%% Eq. (12): 2*M-point IDFT
n = 16; % order of system
% q = 8*n; % Sentence after Eq. (7)
% r = 8*n;

dum_h_hat_i = zeros(p,m,2*M);
for i = 0:q+r-1 % i = 0:2*M-1 This is the correct form but we do not use it
    for k = 0:M-1
        dum_h_hat_i(:,:,i+1) = dum_h_hat_i(:,:,i+1) + (G(:,:,k+1)).*exp(1i*(i-1)*wk(k+1,1)); % exp(1i*(i-1)*wk(k+1,1))  exp(1i*(i-1)*2*pi*k/(2*M))
    end
%     disp(i)
end
h_hat_i = (1/(2*M)).*(dum_h_hat_i);
% h_hat_i = (1/(2*M)).*real(dum_h_hat_i);


% Eq. (13): Assign q and r and construct Hankel matrix H_hat
k = 0;
% H_hat = zeros(q*p,r*m);
for i = 1:q
    for j = 1:r
        H_hat((i-1)*p+1:i*p,(j-1)*m+1:j*m) = squeeze(h_hat_i(:,:,k+1)); % 'k+1' instead of 'k' because i ~= 0 in Eq. (13)
        k = k+1;
    end
    k = i+1; % WRONG: This 'i+2' instead of 'i+1' is because of what I have written in Eq. (12)
end

% SVD of Hankel matrix 
% [U_hat,~,~] = svd(H_hat);
[U_hat,~,~] = svd([real(H_hat) imag(H_hat)]);


% Eq. (14): System order decomposition
U_hat_s = U_hat(:,1:n);

% Eqs. (17)-(19)
J1 = [eye((q-1)*p) zeros((q-1)*p,p)];
J2 = [zeros((q-1)*p,p) eye((q-1)*p)];
J3 = [eye(p) zeros(p,(q-1)*p)];

% Eqs (15) and (16): Finding A_hat and C_hat
% X_plus = @(X) ((X'*X)^(-1))*X'; % Moore-Penrose pseudoinverse of full-column rand matrix X
% A_hat = (X_plus(J1*U_hat_s))*J2*U_hat_s;

% Help pinv
A_hat = (pinv(J1*U_hat_s))*J2*U_hat_s;
C_hat = J3*U_hat_s;

%% Eq. (20): Least Square Problem for detemining B_hat & D_hat 
II = repmat(eye(p),M,1);
S = zeros(M*p,n);
for i = 1:M
    S((i-1)*p+1:i*p,:) = C_hat/((exp(1i*wk(i,1)))*eye(size(A_hat))-A_hat);
end
AA = [S II];
AAA = kron(eye(m),AA);
AAAA = [real(AAA);imag(AAA)];

k = 1;
BBB = zeros(p*M*m,1);
for i = 1:m
    for j = 1:M
        BBB((k-1)*p+1:k*p,1) = squeeze(G(:,i,j));
        k = k+1;
    end
end

BBBB = [real(BBB);imag(BBB)];

% teta = ((AAAA'*AAAA)^-1)*AAAA'*BBBB;
teta = pinv(AAAA)*BBBB;

B_hat = zeros(n,m);
D_hat = zeros(p,m);
for i = 1:m
    teta_dum = teta((i-1)*(p+n)+1:i*(p+n),1);
    B_hat(:,i) = teta_dum(1:n,:);
    D_hat(:,i) = teta_dum(n+1:end,:);
end

%%
epsilon = zeros(p,m,M);
normed_eps = zeros(1,M);
for i = 1:M
    epsilon(:,:,i) = squeeze(G(:,:,i))-D_hat-(C_hat/((exp(1i*wk(i,1)))*eye(size(A_hat))-A_hat))*B_hat; % estimation error
    normed_eps(1,i) = norm(squeeze(epsilon(:,:,i)));
end
ID_error = sum(normed_eps);
