function ID_error = opt_func(inp)

% Fast Fourrier Transformation parameters
% Baseband

m = 1;                 % number of input
p = 1;                 % number of output 
n = 16; % Upper bound of system order

% q = 37; % Sentence after Eq. (7)
% alfa = 60;
q = round(inp(1));
alfa = round(inp(2));

if (n > q) || (n > alfa)
    error('Wrong! q > n OR alfa >= n is not satisfied!!! See step 1 pp. 5')
end


%% Step 1:
load('iddata.mat')
freq_vec = act1_1(:,2);
ww = freq_vec*2*pi;
% omega = freq_vec*2*pi;
T = 5/max(ww);
select_w = 1;
omega = 2*atan(T*ww(select_w:end)/2);

iddata_act1 = (act1_1(:,3)+1i*act1_1(:,4))';
iddata_act2 = (act2_1(:,3)+1i*act2_1(:,4))';
iddata_shaker = (shaker_1(:,3)+1i*shaker_1(:,4))';

% Up to freq(index) Hz 
index = 6401;
freq_vec = freq_vec(1:index,:);
iddata_act1 = iddata_act1(:,1:index);
iddata_act2 = iddata_act2(:,1:index);
iddata_shaker = iddata_shaker(:,1:index);

[N,~] = size(freq_vec); % M-> number of frequency samples

G = zeros(p,m,N);
if ((m == 1)&&(p == 1))
    G(1,1,:) = iddata_act1(1,:);
elseif ((m == 2)&&(p == 1))
    for i = 1:N
        G(1,1,i) = iddata_act1(1,i);
        G(1,2,i) = iddata_shaker(1,i);
    end
else
    for i = 1:N
    G(:,1,i) = iddata_act1(:,i);
    G(:,2,i) = iddata_act2(:,i);
    G(:,3,i) = iddata_shaker(:,i);
    end
end

%% Step 2:
% Create frequency function in Eq. (10)
W_q_N_p_dum = zeros(q,N);
for i = 1:N
   W_q_N_p_dum(:,i) = W_omega(omega(i,1),q);
end
W_q_N_p = kron(W_q_N_p_dum,eye(p,p));

% Create input vector as in Eq. (3)
U = ones(1,m);

% Eq. (4): Create output vector
Y = zeros(p,N);
for i = 1:N
    Y(:,i) = squeeze(G(:,:,i))*U;
end

% Create Y_q_N based on Eq. (11)
Y_q_N_dum = zeros(p*N,N);
for i = 1:N
    Y_q_N_dum((i-1)*p+1:i*p,i) = Y(:,i);
end
Y_q_N = 1/sqrt(N)*W_q_N_p*Y_q_N_dum;

% Eq. (29)
W_q_alfa_N_p_dum = zeros(q+alfa,N);
for i = 1:N
   W_q_alfa_N_p_dum(:,i) = W_omega(omega(i,1),q+alfa);
end
W_q_alfa_N_p = kron(W_q_alfa_N_p_dum,eye(p,p));
U_q_alfa_N_dum = zeros(m*N,N);
for i = 1:N
    U_q_alfa_N_dum((i-1)*m+1:i*m,i) = U;
end
U_q_alfa_N = 1/sqrt(N)*W_q_alfa_N_p*U_q_alfa_N_dum;

%% Step 3: Do the QR factorization
% Last equation in pp.5 
[~,R_dum] = qr([real(U_q_alfa_N) imag(U_q_alfa_N);real(Y_q_N) imag(Y_q_N)]',0); % !!!!!!!! Wrong
R_dum = R_dum';

R_32 = R_dum(m*q+alfa*m+1:end,m*q+1:m*q+alfa*m);


%% Step 4:
[Un,~,~] = svd([real(R_32) imag(R_32)]);

Us = Un(:,1:n);

% Eq. (25)
J1 = [eye((q-1)*p) zeros((q-1)*p,p)];
J2 = [zeros((q-1)*p,p) eye((q-1)*p)];
J3 = [eye(p) zeros(p,(q-1)*p)];

X_plus = @(X) ((X'*X)^(-1))*X'; % Moore-Penrose pseudoinverse of full-column rand matrix X

%% Step 5:
% Eq. (22)
A_hat = (X_plus(J1*Us))*J2*Us;
% Help pinv
% A_hat = (pinv(J1*U_hat_s))*J2*U_hat_s;

% Eq. (23)
C_hat = J3*Us;

%% Step 6:
% Eq. (15): Least Square Problem for detemining B_hat & D_hat 
II = repmat(eye(p),N,1);
S = zeros(N*p,n);
for i = 1:N
    S((i-1)*p+1:i*p,:) = C_hat/((exp(1i*omega(i,1)))*eye(size(A_hat))-A_hat);
end
AA = [S II];
AAA = kron(eye(m),AA);
AAAA = [real(AAA);imag(AAA)];

k = 1;
BBB = zeros(p*N*m,1);
for i = 1:m
    for j = 1:N
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
% estimated transfer function 
epsilon = zeros(p,m,N);
normed_eps = zeros(1,N);
for i = 1:N
    epsilon(:,:,i) = squeeze(G(:,:,i))-D_hat-(C_hat/((exp(1i*omega(i,1)))*eye(size(A_hat))-A_hat))*B_hat; % estimation error
    normed_eps(1,i) = norm(squeeze(epsilon(:,:,i)));
end
ID_error = sum(normed_eps);

