% %% ===========================================================================================
% % Multi-input Multi-output code based on Algorithm 1 in MCKELVEY et al. 1996 
% % in IEEE TAC. The Script is written by ATTA OVEISI at Ruhr-Uni. Bochum in 
% % Germany. Use of this code without permission is not allowed! Contact: atta.oveisi@rub.de
% % ============================================================================================
% 
% %%
% clc
% clear
% close all
% 
% %% Fast Fourrier Transformation parameters
% % Baseband
% fmax = 800;  % frequency span (Hz)
% fs = 2.56*fmax; % Nyquist frequency (Hz)
% 
% dt = 1/fs;  % sampling time (sec)
% 
% m = 1;                 % number of input
% p = 1;                 % number of output 
% 
% %% Input signal generation with window and overlap
% load('iddata.mat')
% freq_vec = act1_1(:,2);
% 
% iddata_act1 = (act1_1(:,3)+1i*act1_1(:,4))';
% iddata_act2 = (act2_1(:,3)+1i*act2_1(:,4))';
% iddata_shaker = (shaker_1(:,3)+1i*shaker_1(:,4))';
% 
% % Up to freq(index) Hz 
% index = 6401;
% freq_vec = freq_vec(1:index,:);
% iddata_act1 = iddata_act1(:,1:index);
% iddata_act2 = iddata_act2(:,1:index);
% iddata_shaker = iddata_shaker(:,1:index);
% 
% [M,~] = size(freq_vec); % M-> number of frequency samples
% 
% G = zeros(p,m,M);
% if ((m == 1)&&(p == 1))
%     G(1,1,:) = iddata_act1(1,:);
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
% %% Uniformly-spaced frequencies
% wk = ((0:M)*pi/M)';
% % wk = 2*atan(dt*(0:M)/M*fmax/2)';
% 
% %% Eq. (11): Extend the transfer function to the full unit circle.
% for k = 1:M-1
%     
%     for i = 1:p
%         for j = 1:m
%             G(i,j,M+k) = conj(G(i,j,M-k)); % if MIMO maybe a transpose is required
%         end
%     end
% end
% 
% %% Eq. (12): 2*M-point IDFT
% dum_h_hat_i = zeros(p,m,2*M);
% for i = 0:2*M-1
%     for k = 0:M-1
%         dum_h_hat_i(:,:,i+1) = dum_h_hat_i(:,:,i+1) + (G(:,:,k+1)).*exp(1i*2*pi*i*k/(2*M));
%     end
% %     disp(i)
% end
% h_hat_i = (1/(2*M)).*dum_h_hat_i;
% 
% save variables

%% Eq. (13): Assign q and r and construct Hankel matrix H_hat
clear
clc

load variables
for i = 16*5:16*10;
    n = 16; % order of system
    q = i; % Sentence after Eq. (7)
    r = q;

    if q+r > 2*M
        error('Wrong! q+r cannot be greater than 2*M')
    end

    if (n > q) || (n > r)
        error('Wrong! q >= n OR r >= n is not satisfied!!!')
    end

    k = 1;
    H_hat = zeros(q*p,r*m);
    for i = 1:q
        for j = 1:r
            H_hat((i-1)*p+1:i*p,(j-1)*m+1:j*m) = squeeze(h_hat_i(:,:,k+1)); % 'k+1' instead of 'k' because i ~= 0 in Eq. (13)
            k = k+1;
        end
        k = i+2; % This 'i+2' instead of 'i+1' is because of what I have written in Eq. (12)
    end

    % SVD of Hankel matrix 
    [U_hat,Sigma_hat,V_hat] = svd(H_hat);
    V_hat_t = V_hat';

    % Eq. (14): System order decomposition
    U_hat_s = U_hat(:,1:n);
    U_hat_o = U_hat(:,n+1:end);

    Sigma_hat_s = Sigma_hat(1:n,1:n);
    Sigma_hat_o = Sigma_hat(n+1:end,n+1:end);

    V_hat_t_s = V_hat_t(1:n,:);
    V_hat_t_o = V_hat_t(n+1:end,:);
    V_hat_s = V_hat_t_s';
    V_hat_o = V_hat_t_o';

    % Eqs. (17)-(19)
    J1 = [eye((q-1)*p) zeros((q-1)*p,p)];
    J2 = [zeros((q-1)*p,p) eye((q-1)*p)];
    J3 = [eye(p) zeros(p,(q-1)*p)];

    % Eqs (15) and (16): Finding A_hat and C_hat
    X_plus = @(X) ((X'*X)^(-1))*X'; % Moore-Penrose pseudoinverse of full-column rand matrix X
    % A_hat = (X_plus(J1*U_hat_s))*J2*U_hat_s;

    % Help pinv
    A_hat = (pinv(J1*U_hat_s))*J2*U_hat_s;
    ss1 = eig(A_hat);
    ss2 = [real(ss1) imag(ss1)];
    figure(1);
    hold on
    if mod(i,3) == 0
        close all
    end
    plot(ss2(:,1),ss2(:,2),'o')
    title(['q and r = ',num2str(q)])
    pause(1)
end

% C_hat = J3*U_hat_s;
% 
% % Eq. (20): Least Square Problem for detemining B_hat & D_hat 
% S = [];
% II = [];
% for i = 1:M+1
%     S = [S ; C_hat*inv((exp(1i*wk(i,1)))*eye(size(A_hat))-A_hat)];
%     II = [II ; eye(p)];
% end
% AA = [S II];
% AAA = [];
% 
% for i = 1:m
%     AAA = blkdiag(AAA,AA);
% end
% 
% AAAA = [real(AAA);imag(AAA)];
% 
% % k = 1;
% % for i = 1:M+1
% %     for j = 1:m
% %         BBB((k-1)*p+1:k*p,1) = squeeze(G(:,j,i));
% %         k = k+1;
% %     end
% % end
% k = 1;
% for i = 1:m
%     for j = 1:M+1
%         BBB((k-1)*p+1:k*p,1) = squeeze(G(:,i,j));
%         k = k+1;
%     end
% end
% 
% BBBB = [real(BBB);imag(BBB)];
% 
% teta = ((AAAA'*AAAA)^-1)*AAAA'*BBBB;
% 
% for i = 1:m
%     teta_dum = teta((i-1)*(p+n)+1:i*(p+n),1);
%     B_hat(:,i) = teta_dum(1:n,:);
%     D_hat(:,i) = teta_dum(n+1:end,:);
% end
% 
% % %%
% % for i = 1:M+1
% %     epsilon(:,:,i) = squeeze(G(:,:,i))-D_hat-(C_hat/((exp(1i*wk(i,1)))*eye(size(A_hat))-A_hat))*B_hat; % estimation error
% %     normed_eps(1,i) = norm(squeeze(epsilon(:,:,i)));
% % end
% % figure(1);
% % plot(wk,normed_eps);
% 
% % estimated transfer function 
% G_hat = ss(A_hat,B_hat,C_hat,D_hat,dt);
% sysc = d2c(G_hat,'zoh');
% figure(2)
% bode(G_hat)
% hold on
% bode(sysc)
% 
% figure(3)
% k = 0;
% for i = 1:m
%     for j = 1:p
%         k = k+1;
%         subplot(m,p,k)
%         plot(freq_vec,20*log10(abs(squeeze(G(j,i,1:M)))));
%         hold on
%         [magG,phaseG,w] = bode(G_hat(j,i),freq_vec*2*pi);
%         for l = 1:size(w)
%             bb(l,1) = squeeze(magG(:,:,l));
%         end
%         plot(w/2/pi,20*log10(bb),'r')
%         xlim([0 fmax])
%     end
% end